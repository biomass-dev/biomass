import numpy as np
from math import fabs, log
from scipy.integrate import simps


def _get_duration(temporal_dynamics):
    """
    Calculation of the duration as the time it takes to decline below 10% of its maximum.

    Parameters
    ----------
    temporal_dynamics: array
        Simulated time course data
    
    Returns
    -------
    duration: int
    
    """
    maximum_value = np.max(temporal_dynamics)
    t_max = np.argmax(temporal_dynamics)
    temporal_dynamics = temporal_dynamics - 0.1 * maximum_value  # 0.1 -> 10 %
    temporal_dynamics[temporal_dynamics > 0.0] = -np.inf
    duration = np.argmax(temporal_dynamics[t_max:]) + t_max

    return duration


def get_signaling_metric(metric, temporal_dynamics):
    """Quantification of cellular response.

    Parameters
    ----------
    metric: str
        'amplitude', 'duration' or 'integral'
    temporal_dynamics: array
        Simulated time course data

    Returns
    -------
    M*: float
        signaling_metric[i, j, k, l]
    
    """
    if metric == 'amplitude':
        return np.max(
            temporal_dynamics
        )
    elif metric == 'duration':
        return _get_duration(
            temporal_dynamics
        )
    elif metric == 'integral':
        return simps(
            temporal_dynamics
        )
    else:
        raise ValueError(
            "Available metrics are: 'amplitude', 'duration', 'integral'"
        )


def dlnyi_dlnxj(signaling_metric, n_file, perturbed_idx,
                observables, conditions, rate, epsilon=1e-9):
    """
    Numerical computation of sensitivities using finite difference approximations
    with 1% changes in the reaction rates or non-zero initial values.

    Parameters
    ----------
    signaling_metric: numpy array
        Signaling metric
    n_file: int
        Number of optimized parameter sets in out/
    perturbed_idx: list
        Indices of rate equations or non-zero initial values
    observables: list
        observables in observable.py
    conditions: list
        Experimental conditions
    rate: float ~ 1
        1.01 for 1% change
    epsilon: float << 1
        If |M - M*| < epsilon, sensitivity_coefficient = 0

    Returns
    -------
    sensitivity_coefficients: numpy array

    """
    sensitivity_coefficients = np.empty(
        (len(n_file), len(perturbed_idx), len(observables), len(conditions))
    )
    for i, _ in enumerate(n_file):
        for j, _ in enumerate(perturbed_idx):
            for k, _ in enumerate(observables):
                for l, _ in enumerate(conditions):
                    if np.isnan(signaling_metric[i, j, k, l]):
                        sensitivity_coefficients[i, j, k, l] = np.nan
                    elif fabs(
                        signaling_metric[i, j, k, l] - signaling_metric[i, -1, k, l]
                    ) < epsilon or (
                        signaling_metric[i, j, k, l] / signaling_metric[i, -1, k, l]
                    ) < 0:
                        sensitivity_coefficients[i, j, k, l] = 0.0
                    else:
                        sensitivity_coefficients[i, j, k, l] = (
                            log(
                                signaling_metric[i, j, k, l] /
                                signaling_metric[i, -1, k, l]
                            ) / log(
                                rate
                            )
                        )

    return sensitivity_coefficients