import os
import re
from types import ModuleType
from typing import List, NamedTuple

import numpy as np

from .graph import NetworkGraph


class OptimizedValues(NamedTuple):
    params: list
    initials: list


class ModelObject(NetworkGraph):
    """
    The BioMASS model object.

    Examples
    --------
    >>> from biomass import Model
    >>> from biomass.models import mapk_cascade
    >>> model = Model(mapk_cascade.__package__).create()
    >>> type(model)
    <class 'biomass.exec_model.ModelObject'>
    >>> print('Parameters:', len(model.parameters))
    Parameters: 22
    >>> print('Species:', len(model.species))
    Species: 8
    >>> print('Observables:', len(model.observables))
    Observables: 2
    """

    def __init__(self, path: str, biomass_model: ModuleType):
        super().__init__(path, biomass_model)

    def get_individual(self, paramset_id: int) -> np.ndarray:
        """
        Get estimated parameter values from optimization results.

        Parameters
        ----------
        paramset_id : int
            Index of parameter set.

        Returns
        -------
        best_individual : numpy.ndarray
            Estimated parameter values.
        """
        best_generation = np.load(
            os.path.join(
                self.path,
                "out",
                f"{paramset_id}",
                "generation.npy",
            )
        )
        best_individual = np.load(
            os.path.join(
                self.path,
                "out",
                f"{paramset_id}",
                f"fit_param{int(best_generation)}.npy",
            )
        )
        return best_individual

    def load_param(self, paramset: int) -> OptimizedValues:
        """
        Load a parameter set from optimization results.

        Parameters
        ----------
        paramset : int
            Index of parameter set.

        Returns
        -------
        optimized_values : OptimizedValues
            Optimized parameter/initial values.
        """
        best_individual = self.get_individual(paramset)
        (x, y0) = self.problem.update(best_individual)
        optimized_values = OptimizedValues(x, y0)
        return optimized_values

    def get_executable(self) -> List[int]:
        """
        Get executable parameter sets from optimization results.
        """
        n_file = []
        try:
            fitparam_files = os.listdir(
                os.path.join(
                    self.path,
                    "out",
                )
            )
            for file in fitparam_files:
                if re.match(r"\d", file):
                    n_file.append(int(file))
            empty_folder = []
            for i, nth_paramset in enumerate(n_file):
                if not os.path.isfile(
                    os.path.join(
                        self.path,
                        "out",
                        f"{nth_paramset:d}",
                        "generation.npy",
                    )
                ):
                    empty_folder.append(i)
            for i in sorted(empty_folder, reverse=True):
                n_file.pop(i)
        except FileNotFoundError as e:
            print(e)
        return n_file

    def get_obj_val(self, indiv_gene: np.ndarray) -> float:
        """
        An objective function to minimize in parameter estimation.

        Parameters
        ----------
        indiv_gene : numpy.ndarray
            Genes, not parameter values.

        Returns
        -------
        obj_val : float
            Objective function value.
        """
        obj_val = self.problem.objective(self.problem.gene2val(indiv_gene))
        return obj_val
