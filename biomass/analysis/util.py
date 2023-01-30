from dataclasses import dataclass, field
from math import isnan, log
from typing import Callable, Dict, Union

import numpy as np
from numba import njit
from scipy.integrate import simpson


@dataclass
class SignalingMetric(object):
    """
    Signaling metric used for sensitivity analysis.
    Attributes
    ----------
    quantification : Dict[str, Callable[[np.ndarray], Union[int, float]]]
        Functions to quantify signaling metrics.
    """

    quantification: Dict[str, Callable[[np.ndarray], Union[int, float]]] = field(
        default_factory=lambda: dict(
            maximum=np.max,
            minimum=np.min,
            integral=simpson,
        ),
        init=False,
    )


@njit(cache=True, fastmath=True)
def dlnyi_dlnxj(
    signaling_metric: np.ndarray,
    num_file_ids: int,
    num_perturbed_ids: int,
    num_observables: int,
    num_conditions: int,
    rate: float,
) -> np.ndarray:
    """
    Numerical computation of sensitivities using finite difference approximations
    with 1% changes in the reaction rates, parameter values or non-zero initial values.
    Parameters
    ----------
    signaling_metric : numpy array
        Signaling metric
    num_file_ids : int
        Number of ptimized parameter sets in out/
    num_perturbed_ids : int
        Number of parameters, rate equations, or non-zero initial conditions to be perturbed.
    num_observables : int
        Number of observables in observable.py
    num_conditions : int
        Number of experimental conditions.
    rate : float ~ 1
        1.01 for 1% change.
    Returns
    -------
    sensitivity_coefficients: numpy array
    """
    EPS = 2.0**-52.0
    sensitivity_coefficients = np.empty(
        (num_file_ids, num_perturbed_ids, num_observables, num_conditions),
    )
    for i in range(num_file_ids):
        for j in range(num_perturbed_ids):
            for k in range(num_observables):
                for l in range(num_conditions):
                    if isnan(signaling_metric[i, j, k, l]):
                        sensitivity_coefficients[i, j, k, l] = np.nan
                    elif (
                        abs(signaling_metric[i, -1, k, l]) < EPS
                        or abs(signaling_metric[i, j, k, l] - signaling_metric[i, -1, k, l]) < EPS
                        or (signaling_metric[i, j, k, l] / signaling_metric[i, -1, k, l]) <= 0
                    ):
                        # 1. Signaling metric before adding perturbation is zero
                        # 2. Absolute change caused by perturbation is too small
                        # 3. Antilogarithm <= 0
                        sensitivity_coefficients[i, j, k, l] = 0.0
                    else:
                        sensitivity_coefficients[i, j, k, l] = log(
                            signaling_metric[i, j, k, l] / signaling_metric[i, -1, k, l]
                        ) / log(rate)

    return sensitivity_coefficients


def remove_nan(sensitivity_matrix: np.ndarray) -> np.ndarray:
    """
    Remove NaN from sensitivity matrix. This function is used for preprocessing of visualizing
    the result of sensitivity analysis through a heatmap.
    Parameters
    ----------
    sensitivity_matrix : numpy.ndarray
        M x N matrix, where M and M are # of parameter sets and # of perturbed objects, respectively.
    """
    EPS = 2.0**-52.0
    nan_idx = []
    for i in range(sensitivity_matrix.shape[0]):
        if np.isnan(sensitivity_matrix[i, :]).any():
            nan_idx.append(i)
        if np.nanmax(np.abs(sensitivity_matrix[i, :])) < EPS:
            sensitivity_matrix[i, :] = np.zeros(sensitivity_matrix.shape[1])
    return np.delete(sensitivity_matrix, nan_idx, axis=0)
