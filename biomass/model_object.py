import os
import re
from types import ModuleType
from typing import List, NamedTuple

import numpy as np


class OptimizedValues(NamedTuple):
    params: list
    initials: list


class ModelObject(object):
    """
    The BioMASS model object.

    Examples
    --------
    >>> from biomass import create_model
    >>> from biomass.models import mapk_cascade
    >>> model = create_model(mapk_cascade.__package__)
    >>> type(model)
    <class 'biomass.model_object.ModelObject'>
    >>> print('Parameters:', len(model.parameters))
    Parameters: 22
    >>> print('Species:', len(model.species))
    Species: 8
    >>> print('Observables:', len(model.observables))
    Observables: 2
    """

    def __init__(self, path: str, biomass_model: ModuleType):
        self._path = path
        self._parameters = biomass_model.C.NAMES
        self._species = biomass_model.V.NAMES
        self.pval = biomass_model.param_values
        self.ival = biomass_model.initial_values
        self.problem = biomass_model.OptimizationProblem()
        self.viz = biomass_model.Visualization()
        self.rxn = biomass_model.ReactionNetwork()

    @property
    def path(self) -> str:
        return self._path

    @property
    def parameters(self) -> List[str]:
        return self._parameters

    @property
    def species(self) -> list:
        return self._species

    @property
    def observables(self) -> List[str]:
        duplicate = [
            name for name in set(self.problem.obs_names) if self.problem.obs_names.count(name) > 1
        ]
        if not duplicate:
            return self.problem.obs_names
        else:
            raise NameError(f"Duplicate observables: {', '.join(duplicate)}")

    def get_individual(self, paramset_id: int) -> np.ndarray:
        """
        Get estimated parameter values from optimization results.

        Parameters
        ----------
        paramset_id : int
            Index of parameter set.

        Returns
        -------
        best_individual : ``numpy.ndarray``
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
        optimized_values : :class:`biomass.model_object.OptimizedValues`
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

    def gene2val(self, indiv_gene: np.ndarray) -> np.ndarray:
        """
        Convert gene to actual value. Mainly used in parameter estimation.

        Parameters
        ----------
        indiv_gene : ``numpy.ndarray``
            Individual gene.

        Returns
        -------
        indiv_values : ``numpy.ndarray``
            Corresponding values.
        """
        bounds = self.problem.get_region()
        indiv_values = 10 ** (indiv_gene * (bounds[1, :] - bounds[0, :]) + bounds[0, :])
        return indiv_values

    def get_obj_val(self, indiv_gene: np.ndarray) -> float:
        """
        An objective function to minimize in parameter estimation.

        Parameters
        ----------
        indiv_gene : ``numpy.ndarray``
            Genes, not parameter values.

        Returns
        -------
        obj_val : float
            Objective function value.
        """
        obj_val = self.problem.objective(self.gene2val(indiv_gene))
        return obj_val
