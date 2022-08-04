import sys
from dataclasses import dataclass, field
from math import fabs, isnan, log
from typing import Callable, Dict, List, Union

import numpy as np
from numba import njit, prange
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


@njit(fastmath=True)
def dlnyi_dlnxj(
    signaling_metric: np.ndarray,
    file_ids:np.ndarray,
    perturbed_ids: np.ndarray,
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

    file_ids : numpy array
        Optimized parameter sets in out/

    perturbed_ids : numpy array
        Indices of rate equations or non-zero initial values.

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
    sensitivity_coefficients = np.empty(
        (len(file_ids), len(perturbed_ids), num_observables, num_conditions),
    )
    for i, _ in enumerate(file_ids):
        for j, _ in enumerate(perturbed_ids):
            for k in range(num_observables):
                for l in range(num_conditions):
                    if np.isnan(signaling_metric[i, j, k, l]):
                        sensitivity_coefficients[i, j, k, l] = np.nan
                    elif (
                        np.abs(signaling_metric[i, -1, k, l]) < sys.float_info.epsilon
                        or np.abs(signaling_metric[i, j, k, l] - signaling_metric[i, -1, k, l])
                        < sys.float_info.epsilon
                        or (signaling_metric[i, j, k, l] / signaling_metric[i, -1, k, l]) <= 0
                    ):
                        # 1. Signaling metric before adding perturbation is zero
                        # 2. Absolute change caused by perturbation is too small
                        # 3. Antilogarithm <= 0
                        sensitivity_coefficients[i, j, k, l] = 0.0
                    else:
                        sensitivity_coefficients[i, j, k, l] = np.log(
                            signaling_metric[i, j, k, l] / signaling_metric[i, -1, k, l]
                        ) / np.log(rate)

    return sensitivity_coefficients


@njit(fastmath=True, parallel=True)
def remove_nan(sensitivity_matrix: np.ndarray) -> np.ndarray:
    """
    Remove NaN from sensitivity matrix. This function is used for preprocessing of visualizing
    the result of sensitivity analysis through a heatmap.

    Parameters
    ----------
    sensitivity_matrix : numpy.ndarray
        M x N matrix, where M and M are # of parameter sets and # of perturbed objects, respectively.
    """ 
    nan_idx = []
    for i in prange(sensitivity_matrix.shape[0]):
        if np.isnan(sensitivity_matrix[i, :]).any():
            nan_idx.append(i)
        if np.nanmax(np.abs(sensitivity_matrix[i, :])) < sys.float_info.epsilon:
            sensitivity_matrix[i, :] = np.zeros(sensitivity_matrix.shape[1])

    return np.delete(sensitivity_matrix, nan_idx, axis=0)