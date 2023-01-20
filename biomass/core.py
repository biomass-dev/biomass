"""BioMASS core functions"""
import os
import warnings
from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
from types import ModuleType
from typing import Callable, Dict, Literal, Optional, Union

import numpy as np
from scipy.optimize import differential_evolution

from .analysis import InitialConditionSensitivity, ParameterSensitivity, ReactionSensitivity
from .dynamics import SignalingSystems
from .estimation import Optimizer
from .model_object import ModelObject
from .version import __version__

__all__ = ["Model", "create_model", "optimize", "run_simulation", "run_analysis"]


class BiomassIndexError(Exception):
    """
    Error in specifying model paramters or species.
    """

    pass


@dataclass
class Model(object):
    """
    Class for BioMASS model construction.

    Attributes
    ----------
    pkg_name: str
        Path (dot-sepalated) to a biomass model directory.
        Use ``__package__``.
    """

    pkg_name: str

    def _load_model(self) -> Optional[ModuleType]:
        try:
            biomass_model = import_module(self.pkg_name)
            return biomass_model
        except ImportError:
            p = Path(self.pkg_name.replace(".", os.sep))
            print(f"cannot import '{p.name}' from '{p.parent}'.")
        return None

    @staticmethod
    def _check_indices(model: ModelObject) -> None:
        """
        Check indices used in a BioMASS model.
        'C' and 'V' must be used for parameters and species, respectively.
        """
        files = ["ode.py", "search_param.py", "observable.py", "reaction_network.py"]
        for file in files:
            with open(os.path.join(model.path, file)) as f:
                lines = f.readlines()
            for i, line in enumerate(lines, start=1):
                loc = f"line {i:d} in {file}"
                msg = ", Use '{}' for the index of a {}."
                if "x[V." in line:
                    raise BiomassIndexError(loc + msg.format("C", "parameter"))
                elif any(map(line.__contains__, ("y[C.", "y0[C.", "dydt[C."))):
                    raise BiomassIndexError(loc + msg.format("V", "species"))

    @staticmethod
    def _check_normalization(model: ModelObject) -> None:
        """
        Check normalization options in ``observable.py``.
        """
        if model.problem.normalization:
            for obs_name in model.observables:
                for key in model.problem.normalization[obs_name]:
                    if key not in ["timepoint", "condition"]:
                        raise KeyError(
                            f"'{key}': invalid key. Should be either 'timepoint' or 'condition'."
                        )
                if (
                    isinstance(model.problem.normalization[obs_name]["timepoint"], int)
                    and not model.problem.t[0]
                    <= model.problem.normalization[obs_name]["timepoint"]
                    <= model.problem.t[-1]
                ):
                    raise ValueError("Normalization timepoint must lie within problem.t.")
                if not model.problem.normalization[obs_name]["condition"]:
                    model.problem.normalization[obs_name]["condition"] = model.problem.conditions
                else:
                    for c in model.problem.normalization[obs_name]["condition"]:
                        if c not in model.problem.conditions:
                            raise ValueError(
                                f"Normalization condition '{c}' is not defined in problem.conditions."
                            )

    @staticmethod
    def _check_visualization_options(model: ModelObject) -> None:
        """
        Check visualization options in ``viz.py``.
        """
        err_msg = (
            "{} in viz.py, len({}) must be equal to or greater than "
            f"len(self.model.problem.conditions) (={len(model.problem.conditions)})."
        )
        singleplotting = model.viz.get_single_observable_options()
        multiplotting = model.viz.get_multiple_observables_options()
        sensitivity_options = model.viz.get_sensitivity_options()
        # Visualization options for time-course simulation (single-observable)
        for i, obs_name in enumerate(model.observables):
            if len(singleplotting[i].cmap) < len(model.problem.conditions):
                raise ValueError(
                    err_msg.format(
                        f"self.single_observable_options for {obs_name}",
                        f"self.single_observable_options[{i}].cmap",
                    )
                )
            elif model.problem.experiments[i] is not None and len(singleplotting[i].shape) < len(
                model.problem.conditions
            ):
                raise ValueError(
                    err_msg.format(
                        f"self.single_observable_options for {obs_name}",
                        f"self.single_observable_options[{i}].shape",
                    )
                )
        # Visualization options for time-course simulation (multi-observables)
        if len(multiplotting.cmap) < len(model.problem.conditions):
            raise ValueError(
                err_msg.format(
                    "self.multiple_observables_options",
                    "self.multiple_observables_options.cmap",
                )
            )
        for obs_name in multiplotting.observables:
            if model.problem.experiments[model.observables.index(obs_name)] is not None and len(
                multiplotting.shape
            ) < len(model.problem.conditions):
                raise ValueError(
                    err_msg.format(
                        "self.multiple_observables_options",
                        "self.multiple_observables_options.shape",
                    )
                )
        # Visualization options for sensitivity analysis results
        if len(sensitivity_options.cmap) < len(model.problem.conditions):
            raise ValueError(
                err_msg.format(
                    "self.sensitivity_options",
                    "self.sensitivity_options.cmap",
                )
            )

    def create(self, show_info: bool = False) -> ModelObject:
        """
        Build a biomass model.

        Parameters
        ----------
        show_info : bool (default: :obj:`False`)
            Set to :obj:`True` to print the information related to model size.

        Returns
        -------
        model : :class:`biomass.model_object.ModelObject`
            The BioMASS model object.

        Examples
        --------
        >>> from biomass import Model
        >>> import your_model
        >>> model = Model(your_model.__package__).create()
        """
        model = ModelObject(self.pkg_name.replace(".", os.sep), self._load_model())
        self._check_indices(model)
        self._check_normalization(model)
        self._check_visualization_options(model)
        if show_info:
            model_name = Path(model.path).name
            print(f"{model_name} information\n" + ("-" * len(model_name)) + "------------")
            print(f"{len(model.species):d} species")
            print(
                f"{len(model.parameters):d} parameters, "
                f"of which {len(model.problem.idx_params):d} to be estimated"
            )
        return model


def create_model(pkg_name: str, show_info: bool = False) -> ModelObject:
    """
    Create a BioMASS model.

    Parameters
    ----------
    pkg_name : str
        Path (dot-sepalated) to a biomass model directory.
    show_info : bool (default: :obj:`False`)
        Set to :obj:`True` to print the information related to model size.

    Returns
    -------
    model : :class:`biomass.model_object.ModelObject`
        The BioMASS model object.

    Examples
    --------
    >>> from biomass import create_model
    >>> import your_model
    >>> model = create_model(your_model.__package__)
    """
    model = Model(pkg_name).create(show_info)
    return model


def optimize(
    model: ModelObject,
    x_id: int,
    *,
    disp_here: bool = False,
    overwrite: bool = False,
    optimizer_options: Optional[dict] = None,
) -> None:
    """
    Estimate model parameters from experimental data.

    Parameters
    ----------
    model : ModelObject
        Model for parameter estimation.
    x_id : int
        Index of parameter set to estimate.
    disp_here : bool (default: :obj:`False`)
        Whether to show the evaluated *objective* at every iteration.
    overwrite : bool (default: :obj:`False`)
        If :obj:`True`, the directory (``x_id/``) will be overwritten.
    optimizer_options : dict, optional
        Keyword arguments to pass to ``scipy.optimize.differential_evolution``.
        For details, please refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html.

    Examples
    --------
    >>> from biomass import create_model, optimize
    >>> from biomass.models import copy_to_current
    >>> copy_to_current("Nakakuki_Cell_2010")
    >>> model = create_model("Nakakuki_Cell_2010")
    >>> optimize(model, x_id=1)

    Notes
    -----

    * Set simulation conditions and the corresponding experimental data in ``observable.py``
    * Define an objective function to be minimized (:func:`objective`) in ``problem.py``
    * Set lower/upper bounds of parameters to be estimated in ``search_param.py``

    """
    if optimizer_options is None:
        optimizer_options = {}
    optimizer_options.setdefault("strategy", "best1bin")
    optimizer_options.setdefault("maxiter", 50)
    optimizer_options.setdefault("popsize", 3)
    optimizer_options.setdefault("tol", 1e-4)
    optimizer_options.setdefault("mutation", 0.1)
    optimizer_options.setdefault("recombination", 0.5)
    optimizer_options.setdefault("disp", True)
    optimizer_options.setdefault("polish", False)
    optimizer_options.setdefault("workers", 1)

    if not optimizer_options["disp"]:
        raise ValueError(
            "Set optimizer_options['disp'] to True. "
            "If you don't want to see the evaluated objective function at every iteration, "
            "set the keyword argument `disp_here` to False."
        )
    if optimizer_options["polish"]:
        warnings.warn(
            "Setting optimizer_options['polish'] to False is highly recommended.",
            UserWarning,
        )

    optimizer = Optimizer(model, differential_evolution, x_id, disp_here, overwrite)
    res = optimizer.minimize(
        model.get_obj_val,
        [(0, 1) for _ in range(len(model.problem.bounds))],
        **optimizer_options,
    )
    param_values = model.gene2val(res.x)
    optimizer.import_solution(param_values)


def run_simulation(
    model: ModelObject,
    *,
    viz_type: str = "original",
    show_all: bool = False,
    stdev: bool = False,
) -> None:
    """
    Simulate ODE model with estimated parameter values.

    Parameters
    ----------
    model : ModelObject
        Model for simulation.
    viz_type : str
        * 'average':
            The average of simulation results with parameter sets in ``out/``.
        * 'best':
            The best simulation result in ``out/``, simulation with `best_fit_param`.
        * 'original':
            Simulation with the default parameters and initial values defined in ``set_model.py``.
        * 'n(=1,2,...)':
            Use the parameter set in ``out/n/``.
        * 'experiment'
            Draw the experimental data written in ``observable.py`` without simulation results.
    show_all : bool (default: :obj:`False`)
        Whether to show all simulation results.
    stdev : bool (default: :obj:`False`)
        If :obj:`True`, the standard deviation of simulated values will be shown
        (only available for 'average' visualization type).

    Examples
    --------
    >>> from biomass import create_model, run_simulation
    >>> from biomass.models import copy_to_current
    >>> copy_to_current("Nakakuki_Cell_2010")
    >>> model = create_model("Nakakuki_Cell_2010")
    >>> run_simulation(
    ...     model,
    ...     viz_type='average',
    ...     show_all=False,
    ...     stdev=True,
    ... )

    """
    if viz_type not in ["best", "average", "original", "experiment"] and not viz_type.isdecimal():
        raise ValueError(
            "Available viz_type are: 'best','average','original','experiment','n(=1, 2, ...)'"
        )

    SignalingSystems(model).simulate_all(
        viz_type=viz_type,
        show_all=show_all,
        stdev=stdev,
    )


def run_analysis(
    model: ModelObject,
    *,
    target: Literal["reaction", "parameter", "initial_condition"],
    metric: str = "integral",
    create_metrics: Optional[Dict[str, Callable[[np.ndarray], Union[int, float]]]] = None,
    style: Literal["barplot", "heatmap"] = "barplot",
    options: Optional[dict] = None,
) -> None:
    """
    Employ sensitivity analysis to identify critical parameters, species or
    reactions in the complex biological network.

    The sensitivity S(y,x) was calculated according to the following equation:
    S(y,x) = d ln(yi) / d ln (xj), where yi is the signaling metric and xj is
    each nonzero species, parameter value or reaction rate.

    Parameters
    ----------
    model : ModelObject
        Model for sensitivity analysis.
    target : Literal["reaction", "parameter", "initial_condition"]
        Where to add a small perturbation to calculate sensitivity coefficients.
    metric : str (default: 'integral')
        A word to specify the signaling metric.
    create_metrics : Dict[str, Callable[[np.ndarray], Union[int, float]]], optional
        Create user-defined signaling metrics.
    style :  Literal["barplot", "heatmap"] (default: 'barplot')
        * 'barplot'
        * 'heatmap'
    options : dict, optional
        * show_indices : bool (default: :obj:`True`)
            (``target`` == 'reaction') Set to :obj:`True` to put reaction index on each bar.

        * excluded_params : list of strings
            (``target`` == 'parameter') List of parameters which are not used for analysis.

        * excluded_initials : list of strings
            (``target`` == 'initial_condition') List of species which are not used for analysis.

        * overwrite : bool (default: :obj:`True`)
            If :obj:`True`, the sensitivity_coefficients/{target}/{metric}.npy file will be overwritten.

    Examples
    --------
    >>> from biomass import create_model, run_analysis
    >>> from biomass.models import copy_to_current
    >>> copy_to_current("Nakakuki_Cell_2010")
    >>> model = create_model("Nakakuki_Cell_2010")

    Parameters

    >>> run_analysis(
    ...     model,
    ...     target='parameter',
    ...     options = {
    ...         'excluded_params': [
    ...             'a', 'Vn', 'Vc', 'Ligand', 'EGF', 'HRG', 'no_ligand'
    ...         ]
    ...     }
    ... )

    Initial condition

    >>> run_analysis(model, target='initial_condition')

    Reaction

    >>> run_analysis(model, target='reaction')

    """

    if options is None:
        options = {}
    options.setdefault("show_indices", True)
    options.setdefault("excluded_params", [])
    options.setdefault("excluded_initials", [])
    options.setdefault("overwrite", True)

    if target == "reaction":
        ReactionSensitivity(model, create_metrics).analyze(
            metric=metric,
            style=style,
            options=options,
        )
    elif target == "parameter":
        ParameterSensitivity(model, create_metrics).analyze(
            metric=metric,
            style=style,
            options=options,
        )
    elif target == "initial_condition":
        InitialConditionSensitivity(model, create_metrics).analyze(
            metric=metric,
            style=style,
            options=options,
        )
    else:
        raise ValueError(
            "Available targets are: '{}".format(
                "', '".join(["reaction", "parameter", "initial_condition"]) + "'."
            )
        )
