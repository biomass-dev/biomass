"""BioMASS core functions"""
import multiprocessing
import warnings

from biomass.dynamics import SignalingSystems
from biomass.ga import GeneticAlgorithmInit, GeneticAlgorithmContinue
from biomass.analysis import (ReactionSensitivity,
                              InitialConditionSensitivity,
                              ParameterSensitivity)


def run_simulation(
        model, 
        viz_type, 
        show_all=False, 
        stdev=False
):
    """
    Simulate ODE model with estimated parameter values.

        Parameters
        ----------
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
        >>> from biomass import run_simulation
        >>> run_simulation(
                Nakakuki_Cell_2010,
                viz_type='average',
                show_all=False, 
                stdev=True
            )
            
        """
    warnings.filterwarnings('ignore')
    if not viz_type in ['best', 'average', 'original', 'experiment'] \
            and not viz_type.isdecimal():
        raise ValueError(
            "Avairable viz_type are: " \
            "'best','average','original','experiment','n(=1, 2, ...)'"
        )
    SignalingSystems(model).simulate_all(
        viz_type=viz_type, show_all=show_all, stdev=stdev
    )


def optimize(
        model, 
        *args, 
        max_generation=10000, 
        allowable_error=0.0
):
    """ 
    Run GA for parameter estimation.

    Paremters
    ---------
    model : module
        Model for parameter estimation
    
    max_generation : int
        Stop if Generation > max_generation
    
    allowable_error : float
        Stop if Best Fitness <= allowable_error
    
    Example
    -------
    >>> from biomass import optimize
    >>> optimize(
            Nakakuki_Cell_2010, 1, 10, max_generation=10000, allowable_error=0.5
        )

    """
    warnings.filterwarnings('ignore')
    ga_init = GeneticAlgorithmInit(
        model,
        max_generation=max_generation,
        allowable_error=allowable_error
    )
    if len(args) == 1:
        ga_init.run(int(args[0]))
    elif len(args) == 2:
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        p = multiprocessing.Pool(processes=n_proc)
        p.map(ga_init.run, range(int(args[0]), int(args[1]) + 1))
        p.close()
    else:
        raise ValueError('too many values to unpack (expected 2)')


def optimize_continue(
        model, 
        *args, 
        max_generation=10000, 
        allowable_error=0.0,
        p0_bounds=[0.1, 10.]
):
    """ 
    Continue running GA from where you stopped in the last parameter search.

    Paremters
    ---------
    model : module
        Model for parameter estimation
    
    max_generation : int
        Stop if Generation > max_generation
    
    allowable_error : float
        Stop if Best Fitness <= allowable_error
    
    p0_bounds : list
        Generate initial population using best parameter values in the last
        parameter search.
            - lower_bound = po_bounds[0] * best_parameter_value
            - upper_bound = p0_bounds[1] * best_parameter_value

    Example
    -------
    >>> from biomass import optimize_continue
    >>> optimize_continue(
            Nakakuki_Cell_2010, 1, 10, max_generation=20000, allowable_error=0.5
        )

    """
    warnings.filterwarnings('ignore')
    ga_continue = GeneticAlgorithmContinue(
        model,
        max_generation=max_generation,
        allowable_error=allowable_error,
        p0_bounds=p0_bounds
    )
    if len(args) == 1:
        ga_continue.run(int(args[0]))
    elif len(args) == 2:
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        p = multiprocessing.Pool(processes=n_proc)
        p.map(ga_continue.run, range(int(args[0]), int(args[1]) + 1))
        p.close()
    else:
        raise ValueError('too many values to unpack (expected 2)')


def run_analysis(
        model,
        target,
        metric='integral',
        style='barplot',
        excluded_params=[]
):
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
    warnings.filterwarnings('ignore')
    if target == 'reaction':
        ReactionSensitivity(model).analyze(metric=metric, style=style)
    elif target == 'initial_condition':
        InitialConditionSensitivity(model).analyze(metric=metric, style=style)
    elif target == 'parameter':
        ParameterSensitivity(
            model, excluded_params
        ).analyze(metric=metric, style=style)
    else:
        raise ValueError(
            "Available targets are: 'reaction', 'initial_condition' , 'parameter'"
        )

