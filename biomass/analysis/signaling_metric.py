import numpy as np
from math import fabs, log
from scipy.integrate import simps
from typing import List, Union, Optional


def _get_duration(
    timecourse: np.ndarray,
    below_threshold: float,
) -> int:
    """
    Calculation of the duration as the time it takes to decline below the threshold

    Parameters
    ----------
    timecourse : array
        Simulated time course data

    below_threshold : float (from 0.0 to 1.0)
        0.1 for 10% of its maximum

    Returns
    -------
    duration : int

    """
    if not 0.0 < below_threshold < 1.0:
        raise ValueError("below_threshold must lie within (0.0, 1.0).")
    maximum_value = np.max(timecourse)
    t_max = np.argmax(timecourse)
    timecourse = timecourse - below_threshold * maximum_value
    timecourse[timecourse > 0.0] = -np.inf
    duration = np.argmax(timecourse[t_max:]) + t_max

    return duration


def get_signaling_metric(
    metric: str,
    timecourse: np.ndarray,
    options: dict,
) -> Optional[Union[int, float]]:
    """Quantification of cellular response.

    Parameters
    ----------
    metric : str
        'maximum', 'minimum', 'duration' or 'integral'

    timecourse : array
        Simulated time course data

    options: dict
        Options for detailed setting of signaling metric.

    Returns
    -------
    M* : int or float
        signaling_metric[i, j, k, l]

    """
    available_metrics = ["maximum", "minimum", "argmax", "argmin", "timepoint", "duration", "integral"]
    if metric not in available_metrics:
        raise ValueError(f"Available metrics are: {', '.join(available_metrics)}")
    elif metric == "maximum":
        return np.max(timecourse)
    elif metric == "minimum":
        return np.min(timecourse)
    elif metric == "argmax":
        return np.argmax(timecourse)
    elif metric == "argmin":
        return np.argmin(timecourse)
    elif metric == "timepoint":
        return timecourse[options["timepoint"]]
    elif metric == "duration":
        return _get_duration(timecourse, options["duration"])
    elif metric == "integral":
        return simps(timecourse)


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
    with 1% changes in the reaction rates, parameter values or non-zero initial values.

    Parameters
    ----------
    signaling_metric : numpy array
        Signaling metric

    n_file : list of integers
        Optimized parameter sets in out/

    perturbed_idx : list of integers
        Indices of rate equations or non-zero initial values.

    observables : list of strings
        observables in observable.py

    conditions : list of strings
        Experimental conditions.

    rate : float ~ 1
        1.01 for 1% change.

    epsilon : float << 1
        If |M - M*| < epsilon, sensitivity_coefficient = 0.

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