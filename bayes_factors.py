"""
Python module for calculating Bayes factors using the method described in 
Dienes 2014 (https://doi.org/10.3389/fpsyg.2014.00781).

(c) Copyright Thomas Akam 2025. Released under the GPL3 License.
"""

import numpy as np
import pylab as plt
from scipy.stats import norm


def bayes_factor(
    data_mean,
    data_SE,
    H1_distribution,
    H1_value=None,
    uniform_max=None,
    uniform_min=None,
    normal_mode=None,
    normal_SD=None,
    half=None,
    H0_value=0,
    plot=False,
    summary=True,
):
    """Calculate a Bayes factor indicating the strength of evidence for an
    alternative hypothesis H1 relative to a null hypothesis H0.

    The data is specified by its mean and standard error of the mean. The sampling
    distribution of the mean (i.e. the distribution of the sample mean around the
    true population mean) is assumed to be normal. Note - this does not require the
    data to be normally distributed, see the Central Limit Theorem.

    The null hypothesis H0 is that the mean takes a specific value, H0_value, which
    is 0 by default.

    The alternative hypothesis H1 is specified as a distribution, which represents a
    prediction about the data mean under the hypothesis.  The type of distribution
    is specified by the H1_distribution parameter and can be 'uniform', 'normal',
    or 'half-normal'.  To set the parameters of the distribution the user can
    provide a single number, H1_value, which should be a reasonable estimate of
    the expected mean of the data under the alternative hypothesis.  This determines
    the parameters of the distribution following recomendations in Dienes 2014 as:
    'uniform':  A uniform distribution between H0_value and H1_value.
    'normal' : A normal distribution with mean=H1_value and SD=(H1_value-H0_value)/2
    'half-normal' : A half-normal distribution with mode=H0 and SD = H1_value-H0_value
     Alternatively the parameters of the hypothesis H1 can also be specified directly,
     see paramter definitions below.  For more information about choosing the alternative
     hypothesis see Dienes 2014 (https://doi.org/10.3389/fpsyg.2014.00781).

     If the plot argument is set to True the data, H1, and H0 distributions are plotted
     for visualisation.  By default the Bayes Factor is printed along with a
     classification of the strength of evidence using the criteria of Lee and Wagenmaker
     2014 (https://doi.org/10.1017/CBO9781139087759), this can be suppressed by setting
     the summary argument to False.

    Parameters:
    data_mean : mean of the observed data.
    data_SE : standard error of the data mean.
    H1_distribution : Alternative hypothesis distribution: 'uniform', 'normal' or 'half-normal'.
    H1_value : Estimate of data mean under alternative hypothesis.
    uniform_max : Maximum value of uniform H1 distribution.
    uniform_min : Minimum value of uniform H1 distribution.
    normal_mode : Mode of the normal or half-normal H1 distribution.
    normal_SD : Standard deviation of normal or half-normal H1 distribution.
    half : whether to use the 'lower' or 'upper' half of the normal distribution.
    HO : value of the mean under null hypothesis, default=0
    plot : if True plot the distributions.
    summary : if True print summary of the result.
    """
    assert H1_distribution in [
        "uniform",
        "normal",
        "half-normal",
    ], "H1_distribution must be 'uniform', 'normal' or 'half-normal"
    # If H1_value is specified use standard H1 distribution parameter choices from Dienes 2014
    if H1_value:
        if H1_distribution == "uniform":
            uniform_max = max(H1_value, H0_value)
            uniform_min = min(H1_value, H0_value)
        elif H1_distribution == "normal":
            normal_mode = H1_value
            normal_SD = np.abs(H1_value - H0_value) / 2
        elif H1_distribution == "half-normal":
            normal_mode = H0_value
            normal_SD = np.abs(H1_value - H0_value)
            half = "upper" if H1_value > H0_value else "lower"
    # Compute data likelihood under null hypotheseis H0
    likelihood_0 = norm.pdf(data_mean - H0_value, scale=data_SE)
    # Compute data likelihood under alternative hypothesis H1
    if H1_distribution == "uniform":
        likelihood_1 = (
            norm.cdf(uniform_max, loc=data_mean, scale=data_SE)
            - norm.cdf(uniform_min, loc=data_mean, scale=data_SE)
        ) / (uniform_max - uniform_min)
    elif H1_distribution == "normal":
        x_min = min(data_mean - 5 * data_SE, normal_mode - 5 * normal_SD)
        x_max = max(data_mean + 5 * data_SE, normal_mode + 5 * normal_SD)
        likelihood_1 = int_normPDF_prod(data_mean, data_SE, normal_mode, normal_SD, x_min, x_max)
    elif H1_distribution == "half-normal":
        if half == "upper":
            x_min = normal_mode
            x_max = normal_mode + 5 * normal_SD
        elif half == "lower":
            x_min = normal_mode - 5 * normal_SD
            x_max = normal_mode
        likelihood_1 = 2 * int_normPDF_prod(
            data_mean, data_SE, normal_mode, normal_SD, x_min, x_max
        )
    # Compute Bayes factor from likelihood ratio.
    bayes_factor = likelihood_1 / likelihood_0
    # Plotting
    if plot:
        plt.figure()
        x_data = np.linspace(data_mean - 5 * data_SE, data_mean + 5 * data_SE, 100)
        plt.fill_between(
            x_data,
            norm.pdf(x_data, loc=data_mean, scale=data_SE),
            label="data",
            alpha=0.5,
        )
        plt.axvline(H0_value, label="H0", c="m")
        if H1_distribution == "uniform":
            height = 1 / (uniform_max - uniform_min)
            plt.fill_between([uniform_min, uniform_max], [height, height], label="H1", alpha=0.5)
        elif H1_distribution == "normal":
            x_H1 = np.linspace(x_min, x_max, 100)
            plt.fill_between(
                x_H1,
                norm.pdf(x_H1, loc=normal_mode, scale=normal_SD),
                label="H1",
                alpha=0.5,
            )
        elif H1_distribution == "half-normal":
            x_H1 = np.linspace(x_min, x_max, 100)
            plt.fill_between(
                x_H1,
                2 * norm.pdf(x_H1, loc=H0_value, scale=normal_SD),
                label="H1",
                alpha=0.5,
            )
        plt.legend()
        plt.ylim(ymin=0)
    # Report evidence strength using criteria of Lee and Wagenmaker 2014.
    if summary:
        supported_H = "H1" if bayes_factor >= 1 else "H0"
        abs_bf = max(bayes_factor, 1 / bayes_factor)
        if abs_bf < 3:
            strength = "Anecdotal"
        elif abs_bf < 10:
            strength = "Moderate"
        elif abs_bf < 30:
            strength = "Strong"
        elif abs_bf < 100:
            strength = "Very strong"
        else:
            strength = "Extreme"
        print(f"Bayes Factor: {bayes_factor :.3g} - {strength} evidence in favour of {supported_H}")
    return bayes_factor


def int_normPDF_prod(m1, sd1, m2, sd2, x_min, x_max):
    """Integrate the product of two normal probability density functions with means
    and standard deviations (m1, sd1) and (m2, sd2) between x_min and x_max. Uses
    the fact that the product of two Gaussian PDFs is a scaled Gaussian PDF."""
    v1 = sd1**2
    v2 = sd2**2
    m = (m1 * v2 + m2 * v1) / (v1 + v2)  # Mean of resulting Gaussian PDF.
    v = (v1 * v2) / (v1 + v2)  # Variance of resulting Gaussian PDF.
    c = norm.pdf(m1 - m2, scale=np.sqrt(v1 + v2))  # Constant to scale resulting PDF.
    return c * norm.cdf(x_max, loc=m, scale=np.sqrt(v)) - c * norm.cdf(
        x_min, loc=m, scale=np.sqrt(v)
    )
