from typing import List, NoReturn, Union

import numpy as np


def initialize_search_param(
    parameters: List[str],
    species: List[str],
    param_values: list,
    initial_values: list,
    estimated_params: List[int],
    estimated_initials: List[int],
) -> Union[np.ndarray, NoReturn]:
    """
    Initialize search_param.

    Parameters
    ----------
    parameters : List[string]
        Names of model parameters.

    species : List[string]
        Names of model species.

    param_values : List[float]
        Parameter values.

    initial_values : List[float]
        Initial values.

    estimated_params : List[int]
        Indices of parameters to be estimated.

    estimated_initials : List[int]
        Indices of initial conditions to be estimated.

    Returns
    -------
    search_param : numpy ndarray
        Numerical values of parameters and/or initial conditions to be estimated.
    """
    if len(estimated_params) != len(set(estimated_params)):
        raise ValueError(
            "Duplicate parameters (C.): {}".format(
                [
                    parameters[idx]
                    for idx in [
                        name for name in set(estimated_params) if estimated_params.count(name) > 1
                    ]
                ]
            )
        )
    elif len(estimated_initials) != len(set(estimated_initials)):
        raise ValueError(
            "Duplicate species (V.): {}".format(
                [
                    species[idx]
                    for idx in [
                        name
                        for name in set(estimated_initials)
                        if estimated_initials.count(name) > 1
                    ]
                ]
            )
        )
    search_param = np.empty(len(estimated_params) + len(estimated_initials))
    for i, j in enumerate(estimated_params):
        search_param[i] = param_values[j]
    for i, j in enumerate(estimated_initials):
        search_param[i + len(estimated_params)] = initial_values[j]

    if np.any(search_param == 0.0):
        message = "search_param must not contain zero."
        for idx in estimated_params:
            if param_values[int(idx)] == 0.0:
                raise ValueError(f'"C.{parameters[int(idx)]}" in idx_params: ' + message)
        for idx in estimated_initials:
            if initial_values[int(idx)] == 0.0:
                raise ValueError(f'"V.{species[int(idx)]}" in idx_initials: ' + message)

    return search_param


def convert_space(
    region: np.ndarray,
    parameters: List[str],
    species: List[str],
    estimated_params: List[int],
    estimated_initials: List[int],
) -> Union[np.ndarray, NoReturn]:
    """
    Convert search space from linear scale to logarithmic scale.

    Parameters
    ----------
    region : numpy ndarray
        Pre-converted search region.

    parameters : List[string]
        Names of model parameters.

    species : List[string]
        Names of model species.

    estimated_params : List[int]
        Indices of parameters to be estimated.

    estimated_initials : List[int]
        Indices of initial conditions to be estimated.

    Returns
    -------
    region : numpy ndarray
        Search bounds for parameters and/or initial conditions to be estimated.
    """
    for i in range(region.shape[1]):
        if np.min(region[:, i]) < 0.0:
            msg = "region[lower_bound, upper_bound] must be positive."
            if i <= len(parameters):
                raise ValueError(f'"C.{parameters[i]}": ' + msg)
            raise ValueError(f'"V.{species[i - len(parameters)]}": ' + msg)
        elif np.min(region[:, i]) == 0 and np.max(region[:, i]) > 0:
            msg = "lower_bound must be larger than 0."
            if i <= len(parameters):
                raise ValueError(f'"C.{parameters[i]}" ' + msg)
            raise ValueError(f'"V.{species[i - len(parameters)]}" ' + msg)
        elif region[1, i] - region[0, i] < 0.0:
            msg = "lower_bound must be smaller than upper_bound."
            if i <= len(parameters):
                raise ValueError(f'"C.{parameters[i]}" : ' + msg)
            raise ValueError(f'"V.{species[i - len(parameters)]}" : ' + msg)
    difference = list(
        set(np.where(np.any(region != 0.0, axis=0))[0])
        ^ set(np.append(estimated_params, [len(parameters) + idx for idx in estimated_initials]))
    )
    if len(difference) > 0:
        msg = "in both search_idx and region"
        for idx in difference:
            if idx <= len(parameters):
                raise ValueError(f'Set "C.{parameters[int(idx)]}" ' + msg)
            raise ValueError(f'Set "V.{species[int(idx - len(parameters))]}" ' + msg)

    search_rgn = region[:, np.any(region != 0.0, axis=0)]

    return np.log10(search_rgn)
