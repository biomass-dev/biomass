"""BioMASS core functions"""
import multiprocessing
import warnings
from typing import Optional, List

from biomass.dynamics import SignalingSystems
from biomass.estimation import GeneticAlgorithmInit, GeneticAlgorithmContinue
from biomass.analysis import (
    ReactionSensitivity,
    InitialConditionSensitivity,
    ParameterSensitivity,
)

__all__ = ["run_simulation", "optimize", "optimize_continue", "run_analysis"]


def run_simulation(model, viz_type: str, show_all: bool = False, stdev: bool = False) -> None:
    """
    Simulate ODE model with estimated parameter values.

        Parameters
        ----------
        model : module
            Model for simulation.

        viz_type : str
            - 'average':
                The average of simulation results with parameter sets in "out/".
            - 'best':
                The best simulation result in "out/", simulation with
                "best_fit_param".
            - 'original':
                Simulation with the default parameters and initial values
                defined in "set_model.py".
            - 'n(=1,2,...)':
                Use the parameter set in "out/n/".
            - 'experiment'
                Draw the experimental data written in observable.py without
                simulation results.

        show_all : bool
            Whether to show all simulation results.

        stdev : bool
            If True, the standard deviation of simulated values will be shown
            (only available for 'average' visualization type).

        Example
        -------
        >>> from biomass.models import Nakakuki_Cell_2010
        >>> from biomass import run_simulation
        >>> run_simulation(
                Nakakuki_Cell_2010,
                viz_type='average',
                show_all=False,
                stdev=True
            )

    """
    warnings.filterwarnings("ignore")
    if not viz_type in ["best", "average", "original", "experiment"] and not viz_type.isdecimal():
        raise ValueError("Available viz_type are: 'best','average','original','experiment','n(=1, 2, ...)'")
    SignalingSystems(model).simulate_all(viz_type=viz_type, show_all=show_all, stdev=stdev)


def optimize(model, *args, options: Optional[dict] = None) -> None:
    """
    Run GA for parameter estimation.

    Paremters
    ---------
    model : module
        Model for parameter estimation.

    options: dict, optional
        popsize : int (default: 5)
            A multiplier for setting the total population size.
            The population has popsize * len(search_param) individuals.

        max_generation : int (default: 10000)
            Stop if Generation > max_generation.

        allowable_error : float (default: 0.0)
            Stop if Best Fitness <= allowable_error.

        local_search_method : str (default: 'mutation')
            Method used in local search. Should be one of
            - 'mutation' : NDM/MGG
            - 'Powell' : Modified Powell method
            - 'DE' : Differential Evolution (strategy: best2bin)

        n_children : int (default: 200)
            (method='mutation') The number of children generated in NDM/MGG.

        workers : int (default: -1 if len(args) == 1 else 1)
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
    >>> optimize(
            Nakakuki_Cell_2010, 1, 10,
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
    options.setdefault("allowable_error", 0.0)
    options.setdefault("local_search_method", "mutation")
    options.setdefault("n_children", 200)
    options.setdefault("workers", -1 if len(args) == 1 else 1)
    options.setdefault("overwrite", False)

    if len(args) == 2 and options["local_search_method"].lower() == "de" and options["workers"] != 1:
        raise AssertionError("daemonic processes are not allowed to have children. Set options['workers'] to 1.")

    warnings.filterwarnings("ignore")
    ga_init = GeneticAlgorithmInit(model, **options)
    if len(args) == 1:
        ga_init.run(int(args[0]))
    elif len(args) == 2:
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        p = multiprocessing.Pool(processes=n_proc)
        p.map(ga_init.run, range(int(args[0]), int(args[1]) + 1))
        p.close()
    else:
        raise ValueError("too many values to unpack (expected 2)")


def optimize_continue(model, *args, options: Optional[dict] = None) -> None:
    """
    Continue running GA from where you stopped in the last parameter search.

    Paremters
    ---------
    model : module
        Model for parameter estimation.

    options: dict, optional
        popsize : int (default: 5)
            A multiplier for setting the total population size.
            The population has popsize * len(search_param) individuals.

        max_generation : int (default: 15000)
            Stop if Generation > max_generation.

        allowable_error : float (default: 0.0)
            Stop if Best Fitness <= allowable_error.

        local_search_method : str (default: 'mutation')
            Method used in local search. Should be one of
            - 'mutation' : NDM/MGG
            - 'Powell' : Modified Powell method
            - 'DE' : Differential Evolution (strategy: best2bin)

        n_children : int (default: 200)
            (method='mutation') The number of children generated in NDM/MGG.

        workers : int (default: -1 if len(args) == 1 else 1)
            (method='DE') The population is subdivided into workers sections and
            evaluated in parallel (uses multiprocessing.Pool). Supply -1 to use
            all available CPU cores. Set workers to 1 when searching multiple
            parameter sets simultaneously.

        p0_bounds : list of float (default: [0.1, 10.0])
            Generate initial population using best parameter values in the last
            parameter search.
                - lower_bound = po_bounds[0] * best_parameter_value
                - upper_bound = p0_bounds[1] * best_parameter_value

    Example
    -------
    >>> from biomass.models import Nakakuki_Cell_2010
    >>> from biomass import optimize_continue
    >>> optimize_continue(
            Nakakuki_Cell_2010, 1, 10,
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
    options.setdefault("allowable_error", 0.0)
    options.setdefault("local_search_method", "mutation")
    options.setdefault("n_children", 200)
    options.setdefault("workers", -1 if len(args) == 1 else 1)
    options.setdefault("p0_bounds", [0.1, 10.0])

    if len(args) == 2 and options["local_search_method"].lower() == "de" and options["workers"] != 1:
        raise AssertionError("daemonic processes are not allowed to have children. Set options['workers'] to 1.")

    warnings.filterwarnings("ignore")
    ga_continue = GeneticAlgorithmContinue(model, **options)
    if len(args) == 1:
        ga_continue.run(int(args[0]))
    elif len(args) == 2:
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        p = multiprocessing.Pool(processes=n_proc)
        p.map(ga_continue.run, range(int(args[0]), int(args[1]) + 1))
        p.close()
    else:
        raise ValueError("too many values to unpack (expected 2)")


def run_analysis(
    model,
    target: str,
    metric: str = "integral",
    style: str = "barplot",
    excluded_params: List[str] = [],
) -> None:
    """
    Perform sensitivity analysis to identify critical parameters, species or
    reactions in the complex biological network.

    The sensitivity S(y,x) was calculated according to the following equation:
    S(y,x) = d ln(yi) / d ln (xj), where yi is the signaling metric and xj is
    each nonzero species, parameter value or reaction rate.

    Paremters
    ---------
    model : module
        Model for sensitivity analysis.

    target : str
        - 'reaction'
        - 'initial_condition'
        - 'parameter'

    metric : str (default: 'integral')
        - 'maximum' : The maximum value.
        - 'minimum' : The minimum value.
        - 'duration' : The time it takes to decline below 10% of its maximum.
        - 'integral' : The integral of concentration over the observation time.

    style : str (default: 'barplot')
        - 'barplot'
        - 'heatmap'

    excluded_params : list of strings
        For parameter sensitivity analysis.

    Example
    -------
    >>> from biomass.models import Nakakuki_Cell_2010
    >>> from biomass import run_analysis
    >>> run_analysis(
            Nakakuki_Cell_2010,
            target='parameter',
            excluded_params=[
                'a', 'Vn', 'Vc', 'Ligand', 'EGF', 'HRG', 'no_ligand'
            ]
        )

    """
    warnings.filterwarnings("ignore")
    if target == "reaction":
        ReactionSensitivity(model).analyze(metric=metric, style=style)
    elif target == "initial_condition":
        InitialConditionSensitivity(model).analyze(metric=metric, style=style)
    elif target == "parameter":
        ParameterSensitivity(model, excluded_params).analyze(metric=metric, style=style)
    else:
        raise ValueError("Available targets are: 'reaction', 'initial_condition' , 'parameter'")
