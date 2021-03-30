"""BioMASS core functions"""
import multiprocessing
import os
from typing import Optional

from .analysis import InitialConditionSensitivity, ParameterSensitivity, ReactionSensitivity
from .dynamics import SignalingSystems
from .estimation import GeneticAlgorithmContinue, GeneticAlgorithmInit
from .exec_model import ModelObject

__all__ = ["optimize", "optimize_continue", "run_simulation", "run_analysis"]


def _check_optional_arguments(
    end: Optional[int],
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
            isinstance(end, int)
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
    start: int,
    end: Optional[int] = None,
    options: Optional[dict] = None,
) -> None:
    """
    Run GA for parameter estimation.

    Paremters
    ---------
    model : ModelObject
        Model for parameter estimation.

    start : int
        Index of parameter set to estimate.

    end : int, optional
        When `end` is specified, parameter sets from `start` to `end` will be estimated.

    options : dict, optional
        popsize : int (default: 5)
            A multiplier for setting the total population size.
            The population has popsize * len(search_param) individuals.

        max_generation : int (default: 10000)
            Stop optimization if Generation > max_generation.

        initial_threshold : float (default: 1e12)
            Threshold on objective function value used to generate initial population.
            Default value is 1e12 (numerically solvable).

        allowable_error : float (default: 0.0)
            Stop optimization if Best Fitness <= allowable_error.

        local_search_method : str (default: 'mutation')
            Method used in local search. Should be one of
            * 'mutation' : NDM/MGG
            * 'Powell' : Modified Powell method
            * 'DE' : Differential Evolution (strategy: best2bin)

        n_children : int (default: 200)
            (method='mutation') The number of children generated in NDM/MGG.

        maxiter : int (default: 10)
            (method='Powell' or 'DE') The maximum number of iterations
            over which the entire population is evolved.

        workers : int (default: -1 if `end` is None else 1)
            (method='DE') The population is subdivided into workers sections and
            evaluated in parallel (uses multiprocessing.Pool). Supply -1 to use
            all available CPU cores. Set workers to 1 when searching multiple
            parameter sets simultaneously.

        overwrite : bool (default: False)
            If True, the out/n folder will be overwritten.

    Example
    -------
    >>> from biomass.models import Nakakuki_Cell_2010
    >>> from biomass import optimize
    >>> model = Nakakuki_Cell_2010.create()
    >>> optimize(
            model=model, start=1, end=10,
            options={
                'max_generation': 10000,
                'allowable_error': 0.5
            }
        )

    """
    if options is None:
        options = {}
    options.setdefault("popsize", 5)
    options.setdefault("max_generation", 10000)
    options.setdefault("initial_threshold", 1e12)
    options.setdefault("allowable_error", 0.0)
    options.setdefault("local_search_method", "mutation")
    options.setdefault("n_children", 200)
    options.setdefault("maxiter", 10)
    options.setdefault("workers", -1 if end is None else 1)
    options.setdefault("overwrite", False)

    _check_optional_arguments(end, options)

    ga_init = GeneticAlgorithmInit(model, **options)
    if end is None:
        ga_init.run(int(start))
    else:
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        with multiprocessing.Pool(processes=n_proc) as p:
            for _ in p.imap_unordered(
                ga_init.run,
                range(int(start), int(end) + 1),
            ):
                pass


def optimize_continue(
    model: ModelObject,
    start: int,
    end: Optional[int] = None,
    options: Optional[dict] = None,
) -> None:
    """
    Continue running GA from where you stopped in the last parameter search.

    Paremters
    ---------
    model : ModelObject
        Model for parameter estimation.

    start : int
        Index of parameter set to estimate.

    end : int, optional
        When `end` is specified, parameter sets from `start` to `end` will be estimated.

    options : dict, optional
        popsize : int (default: 5)
            A multiplier for setting the total population size.
            The population has popsize * len(search_param) individuals.

        max_generation : int (default: 15000)
            Stop optimization if Generation > max_generation.

        initial_threshold : float (default: 1e12)
            Threshold on objective function value used to generate initial population.
            Default value is 1e12 (numerically solvable).

        allowable_error : float (default: 0.0)
            Stop optimization if Best Fitness <= allowable_error.

        local_search_method : str (default: 'mutation')
            Method used in local search. Should be one of
            * 'mutation' : NDM/MGG
            * 'Powell' : Modified Powell method
            * 'DE' : Differential Evolution (strategy: best2bin)

        n_children : int (default: 200)
            (method='mutation') The number of children generated in NDM/MGG.

        maxiter : int (default: 10)
            (method='Powell' or 'DE') The maximum number of iterations
            over which the entire population is evolved.

        workers : int (default: -1 if `end` is None else 1)
            (method='DE') The population is subdivided into workers sections and
            evaluated in parallel (uses multiprocessing.Pool). Supply -1 to use
            all available CPU cores. Set workers to 1 when searching multiple
            parameter sets simultaneously.

        p0_bounds : list of floats (default: [0.1, 10.0])
            Generate initial population using best parameter values in the last
            parameter search.
                - lower_bound = po_bounds[0] * best_parameter_value
                - upper_bound = p0_bounds[1] * best_parameter_value

    Example
    -------
    >>> from biomass.models import Nakakuki_Cell_2010
    >>> from biomass import optimize_continue
    >>> model = Nakakuki_Cell_2010.create()
    >>> optimize_continue(
            model=model, start=1, end=10,
            options={
                'max_generation': 20000,
                'allowable_error': 0.5
            }
        )

    """
    if options is None:
        options = {}
    options.setdefault("popsize", 5)
    options.setdefault("max_generation", 15000)
    options.setdefault("initial_threshold", 1e12)
    options.setdefault("allowable_error", 0.0)
    options.setdefault("local_search_method", "mutation")
    options.setdefault("n_children", 200)
    options.setdefault("maxiter", 10)
    options.setdefault("workers", -1 if end is None else 1)
    options.setdefault("p0_bounds", [0.1, 10.0])

    _check_optional_arguments(end, options)

    ga_continue = GeneticAlgorithmContinue(model, **options)
    if end is None:
        ga_continue.run(int(start))
    else:
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        with multiprocessing.Pool(processes=n_proc) as p:
            for _ in p.imap_unordered(
                ga_continue.run,
                range(int(start), int(end) + 1),
            ):
                pass


def run_simulation(
    model: ModelObject,
    viz_type: str,
    show_all: bool = False,
    stdev: bool = False,
    save_format: str = "pdf",
    param_range: Optional[dict] = None,
) -> None:
    """
    Simulate ODE model with estimated parameter values.

    Parameters
    ----------
    model : ModelObject
        Model for simulation.

    viz_type : str
        * 'average':
            The average of simulation results with parameter sets in "out/".
        * 'best':
            The best simulation result in "out/", simulation with
            "best_fit_param".
        * 'original':
            Simulation with the default parameters and initial values
            defined in "set_model.py".
        * 'n(=1,2,...)':
            Use the parameter set in "out/n/".
        * 'experiment'
            Draw the experimental data written in observable.py without
            simulation results.

    show_all : bool
        Whether to show all simulation results.

    stdev : bool
        If True, the standard deviation of simulated values will be shown
        (only available for 'average' visualization type).

    save_format : str (default: "pdf")
        Either "png" or "pdf", indicating whether to save figures
        as png or pdf format.

    param_range : dict, optional
        orientation : str (default: 'portrait')
            Either 'portrait' or 'landscape'.

        distribution : str (default: 'boxenplot')
            Either 'boxplot' or 'boxenplot'.

        scatter : bool (default: False)
            If True, draw a stripplot.

    Example
    -------
    >>> from biomass.models import Nakakuki_Cell_2010
    >>> from biomass import run_simulation
    >>> model = Nakakuki_Cell_2010.create()
    >>> run_simulation(
            model,
            viz_type='average',
            show_all=False,
            stdev=True,
            save_format="png",
        )

    """
    if viz_type not in ["best", "average", "original", "experiment"] and not viz_type.isdecimal():
        raise ValueError(
            "Available viz_type are: 'best','average','original','experiment','n(=1, 2, ...)'"
        )

    if save_format not in ["pdf", "png"]:
        raise ValueError("save_format must be either 'pdf' or 'png'.")

    if param_range is None:
        param_range = {}
    param_range.setdefault("orientation", "portrait")
    param_range.setdefault("distribution", "boxenplot")
    param_range.setdefault("scatter", False)

    if param_range["orientation"] not in ["portrait", "landscape"]:
        raise ValueError("Available param_range['orientation'] are: 'portrait' or 'landscape'.")
    if param_range["distribution"] not in ["boxplot", "boxenplot"]:
        raise ValueError("Available param_range['distribution'] are: 'boxplot' or 'boxenplot'.")
    if not isinstance(param_range["scatter"], bool):
        raise TypeError("param_range['scatter'] must be a boolean.")

    SignalingSystems(model).simulate_all(
        viz_type=viz_type,
        show_all=show_all,
        stdev=stdev,
        save_format=save_format,
        param_range=param_range,
    )


def run_analysis(
    model: ModelObject,
    target: str,
    metric: str = "integral",
    style: str = "barplot",
    save_format: str = "pdf",
    options: Optional[dict] = None,
) -> None:
    """
    Perform sensitivity analysis to identify critical parameters, species or
    reactions in the complex biological network.

    The sensitivity S(y,x) was calculated according to the following equation:
    S(y,x) = d ln(yi) / d ln (xj), where yi is the signaling metric and xj is
    each nonzero species, parameter value or reaction rate.

    Paremters
    ---------
    model : ModelObject
        Model for sensitivity analysis.

    target : str
        * 'reaction'
        * 'initial_condition'
        * 'parameter'

    metric : str (default: 'integral')
        * 'maximum' : The maximum value.
        * 'minimum' : The minimum value.
        * 'argmax' : The time to reach the maximum value.
        * 'argmin' : The time to reach the minimum value.
        * 'timepoint' : The simulated value at the time point set via options['timepoint'].
        * 'duration' :  The time it takes to decline below the threshold set via options['duration'].
        * 'integral' : The integral of concentration over the observation time.

    style : str (default: 'barplot')
        * 'barplot'
        * 'heatmap'

    save_format : str (default: "pdf")
        Either "png" or "pdf", indicating whether to save figures
        as png or pdf format.

    options : dict, optional
        show_indices : bool (default: True)
            (target == 'reaction') Set to True to put reaction index on each bar.

        excluded_params : list of strings
            (target == 'parameter') List of parameters which are not used for analysis.

        excluded_initials : list of strings
            (target == 'initial_condition') List of species which are not used for analysis.

        timepoint : int (default: model.sim.t[-1])
            (metric=='timepoint') Which timepoint to use.

        duration : float (default: 0.5)
            (metric=='duration') 0.1 for 10% of its maximum.

    Example
    -------
    >>> from biomass.models import Nakakuki_Cell_2010
    >>> from biomass import run_analysis
    >>> model = Nakakuki_Cell_2010.create()

    1. Parameter
    >>> run_analysis(
            model,
            target='parameter',
            options = {
                'excluded_params': [
                    'a', 'Vn', 'Vc', 'Ligand', 'EGF', 'HRG', 'no_ligand'
                ]
            }
        )

    2. Initial condition
    >>> run_analysis(
            model,
            target='initial_condition',
        )

    3. Reaction
    >>> run_analysis(
            model,
            target='reaction',
        )

    """
    if save_format not in ["pdf", "png"]:
        raise ValueError("save_format must be either 'pdf' or 'png'.")

    if options is None:
        options = {}
    options.setdefault("show_indices", True)
    options.setdefault("excluded_params", [])
    options.setdefault("excluded_initials", [])
    options.setdefault("timepoint", model.sim.t[-1])
    options.setdefault("duration", 0.5)

    if not model.sim.t[0] <= options["timepoint"] <= model.sim.t[-1]:
        raise ValueError("options['timepooint'] must lie within sim.t.")
    if not 0.0 < options["duration"] < 1.0:
        raise ValueError("options['duration'] must lie within (0, 1).")

    if target == "reaction":
        ReactionSensitivity(model).analyze(
            metric=metric,
            style=style,
            save_format=save_format,
            options=options,
        )
    elif target == "parameter":
        ParameterSensitivity(model).analyze(
            metric=metric,
            style=style,
            save_format=save_format,
            options=options,
        )
    elif target == "initial_condition":
        InitialConditionSensitivity(model).analyze(
            metric=metric,
            style=style,
            save_format=save_format,
            options=options,
        )
    else:
        here = os.path.abspath(os.path.dirname(__file__))
        files = os.listdir(os.path.join(here, "analysis"))
        raise ValueError(
            "Available targets are: '"
            + "', '".join(
                [
                    available_target
                    for available_target in files
                    if os.path.isdir(
                        os.path.join(
                            here,
                            "analysis",
                            available_target,
                        )
                    )
                    and available_target != "__pycache__"
                ]
            )
            + "'."
        )
