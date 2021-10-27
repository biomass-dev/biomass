"""BioMASS core functions"""
import multiprocessing
import os
from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Optional, Union

try:  # python 3.8+
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

import numpy as np

from .analysis import InitialConditionSensitivity, ParameterSensitivity, ReactionSensitivity
from .dynamics import SignalingSystems
from .estimation import GeneticAlgorithmContinue, GeneticAlgorithmInit
from .exec_model import ModelObject

__all__ = ["Model", "optimize", "optimize_continue", "run_simulation", "run_analysis"]


class BiomassIndexError(Exception):
    """
    Error in specifying model paramters or species.
    """

    pass


@dataclass
class Model(object):
    """
    The BioMASS model object.

    Attributes
    ----------
    pkg_name: str
        Path (dot-sepalated) to a biomass model directory.
        Use ``__package__``.
    """

    pkg_name: str

    def _load_model(self) -> Any:
        try:
            biomass_model = import_module(self.pkg_name)
            return biomass_model
        except ImportError:
            p = Path(self.pkg_name.replace(".", os.sep))
            print(f"cannot import '{p.name}' from '{p.parent}'.")

    @staticmethod
    def _check_indices(model: ModelObject) -> None:
        files = ["set_model.py", "set_search_param.py", "observable.py"]
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
        if model.problem.normalization:
            for obs_name in model.observables:
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

    def create(self, show_info: bool = False) -> ModelObject:
        """
        Build a biomass model.

        Parameters
        ----------
        show_info : bool (default: :obj:`False`)
            Set to :obj:`True` to print the information related to model size.

        Examples
        --------
        >>> from biomass import Model
        >>> import your_model
        >>> model = Model(your_model.__package__).create()
        """
        model = ModelObject(self.pkg_name.replace(".", os.sep), self._load_model())
        self._check_indices(model)
        self._check_normalization(model)
        if show_info:
            model_name = Path(model.path).name
            print(f"{model_name} information\n" + ("-" * len(model_name)) + "------------")
            print(f"{len(model.species):d} species")
            print(
                f"{len(model.parameters):d} parameters, "
                f"of which {len(model.problem.idx_params):d} to be estimated"
            )
        return model


def _check_optional_arguments(
    x_id: Union[int, Iterable[int]],
    options: Optional[dict],
) -> None:
    if options is None:
        pass
    elif isinstance(options, dict):
        if options["local_search_method"].lower() not in ["mutation", "powell", "de"]:
            raise ValueError(
                f"'{options['local_search_method']}': "
                "Invalid local_search_method. Should be one of ['mutation', 'Powell', 'DE']"
            )
        elif (
            not isinstance(x_id, int)
            and options["local_search_method"].lower() == "de"
            and options["workers"] != 1
        ):
            raise AssertionError(
                "daemonic processes are not allowed to have children. Set options['workers'] to 1."
            )
    else:
        raise TypeError("options must be dict or None.")


def optimize(
    model: ModelObject,
    x_id: Union[int, Iterable[int]],
    options: Optional[dict] = None,
) -> None:
    """
    Estimate model parameters from experimental data.

    Parameters
    ----------
    model : ModelObject
        Model for parameter estimation.

    x_id : int or Iterable[int]
        Index (indices) of parameter set to estimate.

    options : dict, optional
        * popsize : int (default: 5)
            A multiplier for setting the total population size.
            The population has ``popsize`` * len(`search_param`) individuals.

        * max_generation : int (default: 10000)
            Stop optimization if Generation > ``max_generation``.

        * initial_threshold : float (default: 1e12)
            Threshold on objective function value used to generate initial population.
            Default value is 1e12 (numerically solvable).

        * allowable_error : float (default: 0.0)
            Stop optimization if Best Fitness <= ``allowable_error``.

        * local_search_method : str (default: 'mutation')
            Method used in local search. Should be one of:

            * 'mutation' : NDM/MGG

            * 'Powell' : Modified Powell method

            * 'DE' : Differential Evolution (strategy: `best2bin`)

        * n_children : int (default: 1000)
            (``local_search_method`` == 'mutation') The number of children generated in NDM/MGG.

        * maxiter : int (default: 10)
            (``local_search_method`` in ['Powell', 'DE']) The maximum number of iterations
            over which the entire population is evolved.

        * workers : int (default: -1 if isinstance(x_id, `int`) else 1)
            (``local_search_method`` == 'DE') The population is subdivided into workers sections and
            evaluated in parallel (uses multiprocessing.Pool). Supply -1 to use
            all available CPU cores. Set workers to 1 when searching multiple
            parameter sets simultaneously.

        * overwrite : bool (default: :obj:`False`)
            If :obj:`True`, the out/n folder will be overwritten.

    Examples
    --------
    >>> from biomass.models import Nakakuki_Cell_2010
    >>> from biomass import Model, optimize
    >>> model = Model(Nakakuki_Cell_2010.__package__).create()
    >>> optimize(
    ...     model, x_id=range(1, 11),
    ...     options={
    ...         'max_generation': 10000,
    ...         'allowable_error': 0.5
    ...     }
    ... )

    """
    if options is None:
        options = {}
    options.setdefault("popsize", 5)
    options.setdefault("max_generation", 10000)
    options.setdefault("initial_threshold", 1e12)
    options.setdefault("allowable_error", 0.0)
    options.setdefault("local_search_method", "mutation")
    options.setdefault("n_children", 1000)
    options.setdefault("maxiter", 10)
    options.setdefault("workers", -1 if isinstance(x_id, int) else 1)
    options.setdefault("overwrite", False)

    _check_optional_arguments(x_id, options)

    ga_init = GeneticAlgorithmInit(model, **options)
    if isinstance(x_id, int):
        ga_init.run(x_id)
    elif all([isinstance(i, int) for i in x_id]):
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        with multiprocessing.Pool(processes=n_proc) as p:
            for _ in p.imap_unordered(ga_init.run, x_id):
                pass
    else:
        raise TypeError("`x_id` type must be either int or Iterable[int].")


def optimize_continue(
    model: ModelObject,
    x_id: Union[int, Iterable[int]],
    options: Optional[dict] = None,
) -> None:
    """
    Continue running optimization from where you stopped in the last parameter search.

    Parameters
    ----------
    model : ModelObject
        Model for parameter estimation.

    x_id : int or Iterable[int]
        Index (indices) of parameter set to estimate.

    options : dict, optional
        * popsize : int (default: 5)
            A multiplier for setting the total population size.
            The population has ``popsize`` * len(`search_param`) individuals.

        * max_generation : int (default: 15000)
            Stop optimization if Generation > ``max_generation``.

        * initial_threshold : float (default: 1e12)
            Threshold on objective function value used to generate initial population.
            Default value is 1e12 (numerically solvable).

        * allowable_error : float (default: 0.0)
            Stop optimization if Best Fitness <= ``allowable_error``.

        * local_search_method : str (default: 'mutation')
            Method used in local search. Should be one of:

            * 'mutation' : NDM/MGG

            * 'Powell' : Modified Powell method

            * 'DE' : Differential Evolution (strategy: `best2bin`)

        * n_children : int (default: 1000)
            (``local_search_method`` == 'mutation') The number of children generated in NDM/MGG.

        * maxiter : int (default: 10)
            (``local_search_method`` in ['Powell', 'DE']) The maximum number of iterations
            over which the entire population is evolved.

        * workers : int (default: -1 if isinstance(x_id, `int`) else 1)
            (``local_search_method`` == 'DE') The population is subdivided into workers sections
            and evaluated in parallel (uses multiprocessing.Pool). Supply -1 to use
            all available CPU cores. Set workers to 1 when searching multiple
            parameter sets simultaneously.

        * p0_bounds : list of floats (default: [0.1, 10.0])
            Generate initial population using best parameter values in the last parameter search.

            * `lower_bound` = po_bounds[0] * `best_parameter_value`
            * `upper_bound` = p0_bounds[1] * `best_parameter_value`

    Examples
    --------
    >>> from biomass.models import Nakakuki_Cell_2010
    >>> from biomass import Model, optimize_continue
    >>> model = Model(Nakakuki_Cell_2010.__package__).create()
    >>> optimize_continue(
    ...     model, x_id=range(1, 11),
    ...     options={
    ...         'max_generation': 20000,
    ...         'allowable_error': 0.5
    ...     }
    ... )

    """
    if options is None:
        options = {}
    options.setdefault("popsize", 5)
    options.setdefault("max_generation", 15000)
    options.setdefault("initial_threshold", 1e12)
    options.setdefault("allowable_error", 0.0)
    options.setdefault("local_search_method", "mutation")
    options.setdefault("n_children", 1000)
    options.setdefault("maxiter", 10)
    options.setdefault("workers", -1 if isinstance(x_id, int) else 1)
    options.setdefault("p0_bounds", [0.1, 10.0])

    _check_optional_arguments(x_id, options)

    ga_continue = GeneticAlgorithmContinue(model, **options)
    if isinstance(x_id, int):
        ga_continue.run(x_id)
    elif all([isinstance(i, int) for i in x_id]):
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        with multiprocessing.Pool(processes=n_proc) as p:
            for _ in p.imap_unordered(ga_continue.run, x_id):
                pass
    else:
        raise TypeError("`x_id` type must be either int or Iterable[int].")


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
    >>> from biomass.models import Nakakuki_Cell_2010
    >>> from biomass import Model, run_simulation
    >>> model = Model(Nakakuki_Cell_2010.__package__).create()
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
    ---------
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

        * overwrite : bool (default: False)
            If :obj:`True`, the sensitivity_coefficients/{target}/{metric}.npy file will be overwritten.

    Examples
    --------
    >>> from biomass.models import Nakakuki_Cell_2010
    >>> from biomass import Model, run_analysis
    >>> model = Model(Nakakuki_Cell_2010.__package__).create()

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
    options.setdefault("overwrite", False)

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
