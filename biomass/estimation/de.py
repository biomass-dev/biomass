import os
import shutil
import sys
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.optimize import OptimizeResult, differential_evolution

from ..exec_model import ExecModel, ModelObject


DIRNAME = "_tmp"


class _Logger(object):
    """
    Duplicate stdout to optimization.log.
    """

    def __init__(self, model_path: str, disp_here: bool):
        self.disp_here = disp_here
        self.terminal = sys.stdout
        self.log_file = open(
            os.path.join(model_path, "out", DIRNAME, "optimization.log"),
            mode="w",
            encoding="utf-8",
        )

    def write(self, message: str):
        if self.disp_here:
            self.terminal.write(message)
        self.log_file.write(message)


@dataclass
class ScipyDifferentialEvolution(ExecModel):
    """
    Use :func:`scipy.optimize.differential_evolution` to estimate kinetic parameters.
    """

    model: ModelObject
    x_id: int

    def __post_init__(self) -> None:
        self.default_stdout = sys.stdout

    def _get_n_iter(self, x_id: int) -> int:

        savedir = os.path.join(self.model.path, "out", f"{x_id:d}")
        n_iter = 0
        with open(
            os.path.join(savedir, "optimization.log"),
            mode="r",
            encoding="utf-8",
        ) as f:
            log_file = f.readlines()
        for message in log_file:
            if len(message.strip()) > 0:
                n_iter += 1
        return n_iter

    def _save_param(self, res: OptimizeResult) -> None:
        """
        Import the solution of the optimization to the model.
        The solution vector `x` will be saved to `path_to_model`/out/`x_id`/.
        Use ``pasmopy.run_simulation`` to visualize the optimization result.
        Parameters
        ----------
        res : OptimizeResult
            The optimization result.
        x_id : int
            Index of the parameter set.
        """
        if os.path.isdir(os.path.join(self.model.path, "out", f"{self.x_id:d}")):
            raise ValueError(
                f"out{os.sep}{self.x_id:d} already exists in {self.model.path}. "
                "Use another parameter id."
            )
        else:
            os.makedirs(os.path.join(self.model.path, "out", f"{self.x_id:d}"))
        savedir = os.path.join(self.model.path, "out", f"{self.x_id:d}")
        shutil.move(os.path.join(self.model.path, "out", DIRNAME, "optimization.log"), savedir)

        param_values = self.model.problem.gene2val(res.x)
        best_fitness: float = self.get_obj_val(res.x)
        n_iter = self._get_n_iter(self.x_id)
        np.save(os.path.join(savedir, "best_fitness"), best_fitness)
        np.save(os.path.join(savedir, "count_num"), n_iter)
        np.save(os.path.join(savedir, "generation"), n_iter)
        np.save(os.path.join(savedir, f"fit_param{n_iter:d}"), param_values)
        # cleanup
        shutil.rmtree(os.path.join(self.model.path, "out", DIRNAME))

    def minimize(self, optimizer_options: dict, disp_here: bool) -> None:
        """
        Run ``scipy.optimize.differential_evolution``.
        The optimization result will be saved in ``model.path/_tmp/x_id``.

        Parameters
        ----------
        x_id: int
            Index  of parameter set to estimate.
        optimizer_options : dict
            Keyword arguments to pass to ``scipy.optimize.differential_evolution``.
            For details, please refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html.
        disp_here: bool
            Whether to show the evaluated *objective* at every iteration.
        Returns
        -------
        res : OptimizeResult
            The optimization result.
        """
        os.makedirs(os.path.join(self.model.path, "out", DIRNAME), exist_ok=True)

        try:
            sys.stdout = _Logger(self.model.path, disp_here)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = differential_evolution(
                    self.get_obj_val,
                    [(0, 1) for _ in range(len(self.model.problem.bounds))],
                    **optimizer_options,
                )
            self._save_param(res)
        finally:
            sys.stdout = self.default_stdout