import numpy as np
from math import fabs, log
from scipy.integrate import simps
from typing import List, Union


def _get_duration(timecourse: np.ndarray, below_threshold: float = 0.1) -> int:
    """
    Calculation of the duration as the time it takes to decline below the threshold

    Parameters
    ----------
    timecourse : array
        Simulated time course data
    below_threshold : float
        0.1 for 10% of its maximum

    Returns
    -------
    duration : int

    """
    maximum_value = np.max(timecourse)
    t_max = np.argmax(timecourse)
    timecourse = timecourse - below_threshold * maximum_value
    timecourse[timecourse > 0.0] = -np.inf
    duration = np.argmax(timecourse[t_max:]) + t_max

    return duration


def get_signaling_metric(metric: str, timecourse: np.ndarray) -> Union[int, float]:
    """Quantification of cellular response.

    Parameters
    ----------
    metric : str
        'maximum', 'minimum', 'duration' or 'integral'
    timecourse : array
        Simulated time course data

    Returns
    -------
    M* : int or float
        signaling_metric[i, j, k, l]

    """
    if metric == "maximum":
        return np.max(timecourse)
    elif metric == "minimum":
        return np.min(timecourse)
    elif metric == "duration":
        return _get_duration(timecourse)
    elif metric == "integral":
        return simps(timecourse)
    else:
        raise ValueError("Available metrics are: 'maximum', 'minimum', 'duration', 'integral'")


def dlnyi_dlnxj(
    signaling_metric: np.ndarray,
    n_file: List[int],
    perturbed_idx: List[int],
    observables: List[str],
    conditions: List[str],
    rate: float,
    epsilon: float = 1e-9,
) -> np.ndarray:
    """
    Numerical computation of sensitivities using finite difference approximations
    with 1% changes in the reaction rates or non-zero initial values.

    Parameters
    ----------
    signaling_metric : numpy array
        Signaling metric
    n_file : list
        Optimized parameter sets in out/
    perturbed_idx : list
        Indices of rate equations or non-zero initial values
    observables : list
        observables in observable.py
    conditions : list
        Experimental conditions
    rate : float ~ 1
        1.01 for 1% change
    epsilon : float << 1
        If |M - M*| < epsilon, sensitivity_coefficient = 0

    Returns
    -------
    sensitivity_coefficients: numpy array

    """
    sensitivity_coefficients = np.empty((len(n_file), len(perturbed_idx), len(observables), len(conditions)))
    for i, _ in enumerate(n_file):
        for j, _ in enumerate(perturbed_idx):
            for k, _ in enumerate(observables):
                for l, _ in enumerate(conditions):
                    if np.isnan(signaling_metric[i, j, k, l]):
                        sensitivity_coefficients[i, j, k, l] = np.nan
                    elif (
                        fabs(signaling_metric[i, j, k, l] - signaling_metric[i, -1, k, l]) < epsilon
                        or (signaling_metric[i, j, k, l] / signaling_metric[i, -1, k, l]) < 0
                    ):
                        sensitivity_coefficients[i, j, k, l] = 0.0
                    else:
                        sensitivity_coefficients[i, j, k, l] = log(
                            signaling_metric[i, j, k, l] / signaling_metric[i, -1, k, l]
                        ) / log(rate)

    return sensitivity_coefficients