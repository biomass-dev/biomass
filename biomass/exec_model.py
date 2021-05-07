import os
import re
from dataclasses import dataclass
from typing import Any, List, NamedTuple

import numpy as np


class ModelObject(object):
    """
    The BioMASS model object.

    Parameters
    ----------
    path : str
        Path to a biomass model.
    biomass_model : Any
        A package containing biomass model properties.

    Attributes
    ----------
    path : str
        Path to the model.
    parameters : list of strings
        Names of model parameters.
    species : list of strings
        Names of model species.
    obs : list of strings
        Names of model observables.
    pval : Callable
        Numerical values of the parameters.
    ival : Callable
        Initial values.
    obj_func : Callable
        An objective function to be minimized for parameter estimation.
    sim : NumericalSimulation
        Simulation conditions and results.
    exp : ExperimentalData
        Experimental measurements used to estimate model parameters.
    viz : Visualization
        Plotting parameters for customizing figure properties.
    sp : SearchParam
        Model parameters to optimize and search region.
    rxn : ReactionNetwork
        Reaction indices grouped according to biological processes.
    """

    def __init__(
        self,
        path: str,
        biomass_model: Any,
    ):
        self._path = path
        self._parameters = biomass_model.C.NAMES
        self._species = biomass_model.V.NAMES
        self._obs = biomass_model.observables
        self.pval = biomass_model.param_values
        self.ival = biomass_model.initial_values
        self.obj_func = biomass_model.objective
        self.sim = biomass_model.NumericalSimulation()
        self.exp = biomass_model.ExperimentalData()
        self.viz = biomass_model.Visualization()
        self.sp = biomass_model.SearchParam()
        self.rxn = biomass_model.ReactionNetwork()

    @property
    def path(self) -> str:
        return self._path

    @property
    def parameters(self) -> List[str]:
        return self._parameters

    @property
    def species(self) -> List[str]:
        return self._species

    @property
    def obs(self) -> List[str]:
        duplicate = [name for name in set(self._obs) if self._obs.count(name) > 1]
        if not duplicate:
            return self._obs
        else:
            raise NameError(f"Duplicate observables: {', '.join(duplicate)}")


class OptimizedValues(NamedTuple):
    params: list
    initials: list


@dataclass
class ExecModel(object):
    model: ModelObject

    def get_individual(self, paramset: int) -> np.ndarray:
        best_generation = np.load(
            os.path.join(
                self.model.path,
                "out",
                f"{paramset:d}",
                "generation.npy",
            )
        )
        best_individual = np.load(
            os.path.join(
                self.model.path,
                "out",
                f"{paramset:d}",
                f"fit_param{int(best_generation):d}.npy",
            )
        )
        return best_individual

    def load_param(self, paramset: int) -> OptimizedValues:
        best_individual = self.get_individual(paramset)
        (x, y0) = self.model.sp.update(best_individual)
        return OptimizedValues(x, y0)

    def get_executable(self) -> List[int]:
        n_file = []
        try:
            fitparam_files = os.listdir(
                os.path.join(
                    self.model.path,
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
                        self.model.path,
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
