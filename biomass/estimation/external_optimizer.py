import os
import shutil
import sys
import warnings
from dataclasses import dataclass
from typing import Callable, List, Union

import numpy as np

from ..exec_model import ModelObject


class _Logger(object):
    """
    Duplicate stdout to _external/optimization.log.
    """

    def __init__(self, model_path: str):
        self.terminal = sys.stdout
        self.log = open(
            os.path.join(model_path, "out", "_external", "optimization.log"),
            mode="w",
            encoding="utf-8",
        )

    def write(self, message: str):
        self.terminal.write(message)
        self.log.write(message)


@dataclass
class ExternalOptimizer(object):
    """
    Use external optimization method for parameterization of a mechanistic model.

    Attributes
    ----------

    model : ModelObject
        The BioMASS model object.

    optimize : Callable
        The external optimizer, e.g., `scipy.optimize.differential_evolution`.
    """

    model: ModelObject
    optimize: Callable

    def __post_init__(self):
        os.makedirs(os.path.join(self.model.path, "out", "_external"), exist_ok=True)
        self.default_stdout = sys.stdout

    def run(self, *args, **kwargs):
        """
        Execute the external optimizer.
        """
        try:
            sys.stdout = _Logger(self.model.path)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = self.optimize(*args, **kwargs)
        finally:
            sys.stdout = self.default_stdout
        return res

    def import_solution(self, x: Union[np.ndarray, List[float]], x_id: int = 0) -> None:
        """
        Import the solution of the optimization to the model.
        The solution vector `x` will be saved to `path_to_model`/out/`x_id`/.
        Use `biomass.run_simulation` to visualize the optimization result.

        Parameters
        ----------
        x : Union[np.ndarray, List[float]]
            The solution of the optimization.
        x_id : int (default: 0)
            Index of the parameter set.

        Examples
        --------
        >>> from scipy.optimize import differential_evolution
        >>> from biomass import Model
        >>> from biomass.models import Nakakuki_Cell_2010
        >>> model = Model(Nakakuki_Cell_2010.__package__).create()
        >>> optimizer = ExternalOptimizer(model, differential_evolution)
        >>> res = optimizer.run(
        ...     model.problem.objective,
        ...     model.problem.bounds,
        ...     strategy="best2bin",
        ...     maxiter=100,
        ...     tol=1e-4,
        ...     mutation=0.1,
        ...     recombination=0.5,
        ...     disp=True,
        ...     polish=False,
        ...     workers=-1,
        ... )

        differential_evolution step 1: f(x)= 7.05589\n
        differential_evolution step 2: f(x)= 5.59166\n
        differential_evolution step 3: f(x)= 2.80301\n
        ...\n
        differential_evolution step 100: f(x)= 0.538524\n

        >>> from biomass import run_simulation
        >>> optimizer.import_solution(res.x, x_id=0)
        >>> run_simulation(model, viz_type="0")
        """

        if os.path.isdir(os.path.join(self.model.path, "out", f"{x_id:d}")):
            raise ValueError(
                f"out{os.sep}{x_id:d} already exists in {self.model.path}. "
                "Use another parameter id."
            )
        else:
            os.makedirs(os.path.join(self.model.path, "out", f"{x_id:d}"))
        shutil.move(
            os.path.join(self.model.path, "out", "_external", "optimization.log"),
            os.path.join(self.model.path, "out", f"{x_id:d}"),
        )

        best_fitness: float = self.model.problem.objective(x)
        n_iter: int = 0
        with open(
            os.path.join(self.model.path, "out", f"{x_id:d}", "optimization.log"),
            mode="r",
            encoding="utf-8",
        ) as f:
            log_file = f.readlines()
        for message in log_file:
            if len(message.strip()) > 0:
                n_iter += 1
        np.save(
            os.path.join(self.model.path, "out", f"{x_id:d}", "best_fitness"),
            best_fitness,
        )
        np.save(
            os.path.join(self.model.path, "out", f"{x_id:d}", "count_num"),
            n_iter,
        )
        np.save(
            os.path.join(self.model.path, "out", f"{x_id:d}", "generation"),
            n_iter,
        )
        np.save(
            os.path.join(self.model.path, "out", f"{x_id:d}", f"fit_param{n_iter:d}"),
            x,
        )
