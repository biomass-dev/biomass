from typing import List, NoReturn

import numpy as np


def initialize_search_param(
    parameters: List[str],
    species: List[str],
    param_values: list,
    initial_values: list,
    estimated_params: List[int],
    estimated_initials: List[int],
) -> np.ndarray:
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
                ", ".join(
                    [
                        parameters[idx]
                        for idx in [
                            name
                            for name in set(estimated_params)
                            if estimated_params.count(name) > 1
                        ]
                    ]
                )
            )
        )

    if len(estimated_initials) != len(set(estimated_initials)):
        raise ValueError(
            "Duplicate initial conditions (V.): {}".format(
                ", ".join(
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
        )

    for name in [parameters[i] for i in estimated_params]:
        if name.startswith("init_") and name.lstrip("init_") in [
            species[j] for j in estimated_initials
        ]:
            raise NameError(
                f"{name}: The parameter name prefix 'init_' "
                "indicates the initial value of a dynamic variable.",
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
                raise ValueError(f"'C.{parameters[int(idx)]}' in self.idx_params: " + message)
        for idx in estimated_initials:
            if initial_values[int(idx)] == 0.0:
                raise ValueError(f"'V.{species[int(idx)]}' in self.idx_initials: " + message)
    return search_param


def _raise_value_error(
    parameters: List[str],
    species: List[str],
    search_param_index: int,
    message: str,
) -> NoReturn:
    """
    Raise ValueError when user-defined search_rgn is not appropriate.
    """
    if search_param_index <= len(parameters):
        raise ValueError(f"'C.{parameters[search_param_index]}': " + message)
    raise ValueError(f"'V.{species[search_param_index - len(parameters)]}': " + message)


def convert_scale(
    region: np.ndarray,
    parameters: List[str],
    species: List[str],
    estimated_params: List[int],
    estimated_initials: List[int],
) -> np.ndarray:
    """
    Convert search region from linear scale to logarithmic scale.

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
    for idx_search in range(region.shape[1]):
        if np.min(region[:, idx_search]) < 0.0:
            _raise_value_error(
                parameters,
                species,
                idx_search,
                "search_rgn[lower_bound, upper_bound] must be positive.",
            )
        elif np.min(region[:, idx_search]) == 0 and np.max(region[:, idx_search]) > 0:
            _raise_value_error(
                parameters,
                species,
                idx_search,
                "lower_bound must be larger than 0.",
            )
        elif region[0, idx_search] > region[1, idx_search]:
            _raise_value_error(
                parameters,
                species,
                idx_search,
                "lower_bound must be smaller than upper_bound.",
            )
    difference = list(
        set(np.where(np.any(region != 0.0, axis=0))[0])
        ^ set(np.append(estimated_params, [len(parameters) + idx for idx in estimated_initials]))
    )
    if len(difference) > 0:
        for idx_diff in difference:
            _raise_value_error(
                parameters,
                species,
                idx_diff,
                "should be set in both search_idx and search_rgn.",
            )

    search_rgn = region[:, np.any(region != 0.0, axis=0)]

    return np.log10(search_rgn)
