import numpy as np
from scipy.integrate import solve_ivp
from typing import List, Callable, Optional

from .name2idx import C, V
from .set_model import DifferentialEquation


observables = []


class NumericalSimulation(DifferentialEquation):
    """Simulate a model using scipy.integrate.ode

    Attributes
    ----------
    normalization : nested dict
        Keys for each observable
        ------------------------
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If None, the maximum value will be used for normalization.

        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in sim.conditions will be used.

    """

    def __init__(self):
        super().__init__(perturbation={})
        self.normalization = {}

    t = range(101)  # 0, 1, 2, ..., 100

    # Experimental conditions
    conditions = []

    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x: list, y0: list, _perturbation: dict = {}) -> Optional[bool]:
        if _perturbation:
            self.perturbation = _perturbation
        # unperturbed steady state

        for i, condition in enumerate(self.conditions):

            sol = self._solveode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                pass

    @staticmethod
    def _solveode(
        diffeq: Callable,
        y0: List[float],
        t: range,
        f_params: tuple,
        method: str = "BDF",
        options: Optional[dict] = None,
    ):
        """
        Solve a system of ordinary differential equations
        using scipy.integrate.solve_ivp.

        Parameters
        ----------
        diffeq : callable f(t, y, *x)
            Right-hand side of the differential equation.

        y0 : array
            Initial condition on y (can be a vector).

        t : array
            A sequence of time points for which to solve for y.

        f_params : tuple
            Model parameters: tuple(x).

        method : str (default: "BDF")
            Integration method to use.

        options : dict, optional
            Options passed to a chosen solver.

        Returns
        -------
        sol : OdeResult
            Represents the solution of ODE.

        """
        if options is None:
            options = {}
        options.setdefault("rtol", 1e-8)
        options.setdefault("atol", 1e-8)
        try:
            sol = solve_ivp(
                diffeq,
                (t[0], t[-1]),
                y0,
                method=method,
                t_eval=t,
                args=f_params,
                **options,
            )
            return sol if sol.success else None
        except ValueError:
            return None

    def _get_steady_state(
        self,
        diffeq: Callable,
        y0: List[float],
        f_params: tuple,
        eps: float = 1e-6,
    ):
        """
        Find the steady state for the untreated condition.

        Parameters
        ----------
        diffeq : callable f(t, y, *x)
            Right-hand side of the differential equation.

        y0 : array
            Initial condition on y (can be a vector).

        f_params : tuple
            Model parameters.

        eps : float (default: 1e-6)
            Run until a time t for which the maximal absolute value of the
            regularized relative derivative was smaller than eps.

        Returns
        -------
        y0 : array
            Steady state concentrations of all species.

        """
        while True:
            sol = self._solveode(diffeq, y0, range(2), f_params)
            if sol is None or np.max(np.abs((sol.y[:, -1] - y0) / (np.array(y0) + eps))) < eps:
                break
            else:
                y0 = sol.y[:, -1].tolist()

        return [] if sol is None else sol.y[:, -1].tolist()


class ExperimentalData(object):
    """
    Set experimental data.

    Attributes
    ----------
    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """

    def __init__(self):
        self.experiments = [None] * len(observables)
        self.error_bars = [None] * len(observables)

    def set_data(self):
        pass

    @staticmethod
    def get_timepoint(obs_name: str) -> Optional[List[int]]:
        if obs_name in observables:
            return []