from math import fabs, log
import numpy as np


def get_duration(time_course_vector):
    """Calculation of the duration as the time it takes to decline below 10% of its maximum
    """
    maximum_value = np.max(time_course_vector)
    t_max = np.argmax(time_course_vector)
    time_course_vector = time_course_vector - 0.1*maximum_value
    time_course_vector[time_course_vector > 0.0] = -np.inf
    duration = np.argmax(time_course_vector[t_max:]) + t_max

    return duration


def compute_sensitivity_coefficients(signaling_metric, n_file, perturbed_idx,
                                        observables, conditions, rate, metric_idx,
                                        epsilon=1e-9):
    """Numerical computation of sensitivities

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
    metric_idx: int (0 or -1)
        0 -> rate equations
        -1 -> non-zero initial values
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
                        signaling_metric[i, j, k, l] - signaling_metric[i, metric_idx, k, l]
                    ) < epsilon or (
                        signaling_metric[i, j, k, l] / signaling_metric[i, metric_idx, k, l]
                    ) < 0:
                        sensitivity_coefficients[i, j, k, l] = 0.0
                    else:
                        sensitivity_coefficients[i, j, k, l] = (
                            log(
                                signaling_metric[i, j, k, l] /
                                signaling_metric[i, metric_idx, k, l]
                            ) / log(rate)
                        )

    return sensitivity_coefficients