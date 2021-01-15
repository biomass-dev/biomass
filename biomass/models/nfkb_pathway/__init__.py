from typing import List

from .name2idx import C, V
from .set_model import param_values, initial_values
from .observable import observables, NumericalSimulation, ExperimentalData
from .viz import Visualization
from .fitness import objective
from .set_search_param import SearchParam
from .reaction_network import ReactionNetwork


def _check_duplicate(names: List[str], object: str) -> List[str]:
    duplicate = [name for name in set(names) if names.count(name) > 1]
    if not duplicate:
        return names
    else:
        raise NameError(f"Duplicate {object}: {', '.join(duplicate)}")


class Model(object):
    def __init__(self) -> None:
        self.path = __path__[0]
        self.parameters = _check_duplicate(C.NAMES, "parameters")
        self.species = _check_duplicate(V.NAMES, "species")
        self.obs = _check_duplicate(observables, "observables")
        self.pval = param_values
        self.ival = initial_values
        self.obj_func = objective
        self.sim = NumericalSimulation()
        self.exp = ExperimentalData()
        self.viz = Visualization()
        self.sp = SearchParam()
        self.rxn = ReactionNetwork()


def show_properties() -> None:
    print(
        "Model properties\n"
        "----------------\n"
        f"{len(Model().species):d} species\n"
        f"{len(Model().parameters):d} parameters, "
        f"of which {len(Model().sp.idx_params):d} to be estimated"
    )


def create() -> Model:
    model = Model()
    if model.sim.normalization:
        for obs_name in model.obs:
            if (
                isinstance(model.sim.normalization[obs_name]["timepoint"], int)
                and not model.sim.t[0] <= model.sim.normalization[obs_name]["timepoint"] <= model.sim.t[-1]
            ):
                raise ValueError("Normalization timepoint must lie within sim.t.")
            if not model.sim.normalization[obs_name]["condition"]:
                model.sim.normalization[obs_name]["condition"] = model.sim.conditions
            else:
                for c in model.sim.normalization[obs_name]["condition"]:
                    if not c in model.sim.conditions:
                        raise ValueError(f"Normalization condition '{c}' is not defined in sim.conditions.")
    return model
