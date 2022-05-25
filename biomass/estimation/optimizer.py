import multiprocessing
import os
import shutil
import sys
import warnings
from dataclasses import dataclass
from math import isfinite
from typing import Callable, List, Union

import numpy as np
from tqdm import tqdm

from ..exec_model import ModelObject

DIRNAME = "_tmp"


class _Logger(object):
    """
    Duplicate stdout to _tmp/optimization.log.
    """

    def __init__(self, model_path: str, x_id: int, disp_here: bool):
        self.disp_here = disp_here
        self.terminal = sys.stdout
        self.log = open(
            os.path.join(model_path, "out", DIRNAME + str(x_id), "optimization.log"),
            mode="w",
            encoding="utf-8",
        )

    def write(self, message: str):
        if self.disp_here:
            self.terminal.write(message)
        self.log.write(message)


class Optimizer(object):
    """
    Use an external optimization method for parameterization of a mechanistic model.

    Attributes
    ----------
    model : ModelObject
        The BioMASS model object.
    optimize : Callable
        The optimizer, e.g., :func:`scipy.optimize.differential_evolution`.
    x_id : int
        Index of parameter set to estimate.
    disp_here: bool (default: False)
        Whether to show the evaluated *objective* at every iteration.

    Examples
    --------
    >>> from scipy.optimize import differential_evolution
    >>> from biomass import Model
    >>> from biomass.estimation import ExternalOptimizer
    >>> from biomass.models import Nakakuki_Cell_2010
    >>> model = Model(Nakakuki_Cell_2010.__package__).create()
    >>> param_idx = 1
    >>> optimizer = ExternalOptimizer(model, differential_evolution, param_idx)
    >>> def obj_fun(x):
    ...    '''Objective function to be minimized.'''
    ...    return optimizer.get_obj_val(x)
    >>> res = optimizer.minimize(
    ...     obj_fun,
    ...     [(0, 1) for _ in range(len(model.problem.bounds))],
    ...     strategy="best1bin",
    ...     maxiter=50,
    ...     tol=1e-4,
    ...     mutation=0.1,
    ...     recombination=0.5,
    ...     disp=True,
    ...     polish=False,
    ...     workers=-1,
    ... )

    differential_evolution step 1: f(x)= 5.19392\n
    differential_evolution step 2: f(x)= 2.32477\n
    differential_evolution step 3: f(x)= 1.93583\n
    ...\n
    differential_evolution step 50: f(x)= 0.519774\n

    >>> from biomass import run_simulation
    >>> param_values = model.problem.gene2val(res.x)
    >>> optimizer.import_solution(param_values)
    >>> run_simulation(model, viz_type=str(param_idx))
    """

    def __init__(
        self,
        model: ModelObject,
        optimize: Callable,
        x_id: int,
        disp_here: bool = False,
    ):
        self.model = model
        self.optimize = optimize
        self.x_id = x_id
        self.disp_here = disp_here

        self.savedir = os.path.join(self.model.path, "out", f"{self.x_id}")
        if os.path.isdir(self.savedir):
            raise ValueError(
                f"out{os.sep}{self.x_id} already exists in {self.model.path}. "
                "Use another parameter id."
            )
        else:
            os.makedirs(self.savedir)
        os.makedirs(os.path.join(self.model.path, "out", DIRNAME + str(self.x_id)), exist_ok=True)
        self.default_stdout = sys.stdout

    def minimize(self, *args, **kwargs):
        """
        Execute the external optimizer.
        """
        os.makedirs(os.path.join(self.model.path, "out", DIRNAME + str(self.x_id)), exist_ok=True)
        try:
            sys.stdout = _Logger(self.model.path, self.x_id, self.disp_here)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = self.optimize(*args, **kwargs)
            return res
        finally:
            sys.stdout = self.default_stdout

    def _get_n_iter(self) -> int:
        n_iter: int = 0
        path_to_log = os.path.join(self.savedir, "optimization.log")
        with open(path_to_log, mode="r", encoding="utf-8") as f:
            log_file = f.readlines()
        for message in log_file:
            if len(message.strip()) > 0:
                n_iter += 1
        return n_iter

    def import_solution(self, x: Union[np.ndarray, List[float]], cleanup: bool = True) -> None:
        """
        Import the solution of the optimization to the model.
        The solution vector `x` will be saved to `path_to_model`/out/`x_id`/.
        Use ``biomass.run_simulation`` to visualize the optimization result.

        Parameters
        ----------
        x : Union[np.ndarray, List[float]]
            The solution of the optimization.
        cleanup : bool (default: True)
            If True (default), delete the temporary folder after the optimization is finished.
        """
        shutil.move(
            os.path.join(self.model.path, "out", DIRNAME + str(self.x_id), "optimization.log"),
            self.savedir,
        )

        best_fitness: float = self.model.problem.objective(x)
        n_iter = self._get_n_iter()
        np.save(os.path.join(self.savedir, "best_fitness"), best_fitness)
        np.save(os.path.join(self.savedir, "count_num"), n_iter)
        np.save(os.path.join(self.savedir, "generation"), n_iter)
        np.save(os.path.join(self.savedir, f"fit_param{n_iter}"), x)
        if cleanup:
            shutil.rmtree(os.path.join(self.model.path, "out", DIRNAME + str(self.x_id)))


@dataclass
class InitialPopulation(object):
    """
    Attributes
    ----------
    model : ModelObject
        BioMASS model object.
    popsize : int (default: 3)
        A multiplier for setting the total population size.
    threshold : float (default: 1e12)
        Allowable error for generating initial population.
        Default value is 1e12 (numerically solvable).
    """

    model: ModelObject
    popsize: int = 3
    threshold: float = 1e12

    def __post_init__(self) -> None:
        self.n_gene = len(self.model.problem.bounds)
        self.n_population = self.popsize * self.n_gene
        self.initpop_path = os.path.join(self.model.path, "_initpop")
        os.makedirs(self.initpop_path, exist_ok=True)

    def _set_gene_vector(self, pop_id: int) -> None:
        obj_val = np.inf
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            while self.threshold <= obj_val or not isfinite(obj_val):
                gene = np.random.rand(self.n_gene)
                obj_val = self.model.get_obj_val(gene)
            np.save(os.path.join(self.initpop_path, f"population_{pop_id}"), gene)

    def generate(self, n_proc: int = 1, progress: bool = False) -> np.ndarray:
        """
        Return initial population for ``optimizer_options['init'] `` in ``biomass.optimize`` function.

        Parameters
        ----------
        n_proc : int (default: 1)
            Number of processes to use (default: 1). Set to a larger number
            (e.g. the number of CPU cores available) for parallel execution of simulations.
        progress : bool (default: :obj:`False`)
            Whether the progress bar is animating or not.

        Returns
        -------
        population : numpy.ndarray
            Array specifying the initial population.
            The array should have shape (M, len(x)),
            where M is the total population size and len(x) is the number of parameters.

        Examples
        --------
        >>> from biomass import create_model, optimize
        >>> from biomass.estimation import InitialPopulation
        >>> from biomass.models import Nakakuki_Cell_2010
        >>> model = create_model(Nakakuki_Cell_2010.__package__)
        >>> initpop = InitialPopulation(model).generate(n_proc=2, progress=True)
        >>> optimize(model, x_id=1, optimizer_options={"init": initpop})
        """
        p = multiprocessing.Pool(processes=n_proc)
        population = np.empty((self.n_population, self.n_gene))
        try:
            with tqdm(total=self.n_population, disable=not progress) as t:
                for _ in p.imap_unordered(self._set_gene_vector, range(self.n_population)):
                    t.update(1)
            for i in range(self.n_population):
                gene = np.load(os.path.join(self.initpop_path, f"population_{i}.npy"))
                population[i] = gene
            return population
        finally:
            p.close()
            shutil.rmtree(self.initpop_path)
