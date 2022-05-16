import os
import shutil
import sys
import warnings
from typing import Callable, List, Union

import numpy as np

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
