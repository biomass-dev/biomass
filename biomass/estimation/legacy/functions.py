import multiprocessing
import warnings
from typing import Iterable, Optional, Union

from ...exec_model import ModelObject
from ._ga import GeneticAlgorithmContinue, GeneticAlgorithmInit


def _check_optional_arguments(
    x_id: Union[int, Iterable[int]],
    options: Optional[dict],
) -> None:
    warnings.warn(
        "This function will no longer be supported. Use biomass.optimize() for parameter estimation.",
        FutureWarning,
    )

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
