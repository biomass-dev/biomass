import os
from typing import List

from .fitness import objective
from .name2idx import C, V
from .observable import ExperimentalData, NumericalSimulation, observables
from .reaction_network import ReactionNetwork
from .set_model import initial_values, param_values
from .set_search_param import SearchParam
from .viz import Visualization


def _check_duplicate(names: List[str], object: str) -> List[str]:
    duplicate = [name for name in set(names) if names.count(name) > 1]
    if not duplicate:
        return names
    else:
        raise NameError(f"Duplicate {object}: {', '.join(duplicate)}")


class BioMassModel(object):
    def __init__(self) -> None:
        self._path = __path__[0]
        self._parameters = _check_duplicate(C.NAMES, "parameters")
        self._species = _check_duplicate(V.NAMES, "species")
        self._obs = _check_duplicate(observables, "observables")
        self.pval = param_values
        self.ival = initial_values
        self.obj_func = objective
        self.sim = NumericalSimulation()
        self.exp = ExperimentalData()
        self.viz = Visualization()
        self.sp = SearchParam()
        self.rxn = ReactionNetwork()

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
        return self._obs


def show_info() -> None:
    model_name = os.path.basename(os.path.dirname(__file__))
    print(
        f"{model_name} information\n" + ("-" * len(model_name)) + "------------\n"
        f"{len(BioMassModel().species):d} species\n"
        f"{len(BioMassModel().parameters):d} parameters, "
        f"of which {len(BioMassModel().sp.idx_params):d} to be estimated"
    )


def create() -> BioMassModel:
    model = BioMassModel()
    if model.sim.normalization:
        for obs_name in model.obs:
            if (
                isinstance(model.sim.normalization[obs_name]["timepoint"], int)
                and not model.sim.t[0]
                <= model.sim.normalization[obs_name]["timepoint"]
                <= model.sim.t[-1]
            ):
                raise ValueError("Normalization timepoint must lie within sim.t.")
            if not model.sim.normalization[obs_name]["condition"]:
                model.sim.normalization[obs_name]["condition"] = model.sim.conditions
            else:
                for c in model.sim.normalization[obs_name]["condition"]:
                    if c not in model.sim.conditions:
                        raise ValueError(
                            f"Normalization condition '{c}' is not defined in sim.conditions."
                        )
    return model
