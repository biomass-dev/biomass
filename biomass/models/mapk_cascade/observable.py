import numpy as np
from scipy.integrate import solve_ivp
from typing import List, Callable, Optional

from .name2idx import C, V
from .set_model import DifferentialEquation


observables = [
    "biphosphorylated_MAPK",
    "unphosphorylated_MAPK",
]


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

    t = range(150 * 60 + 1)

    # Experimental conditions
    conditions = ["control"]

    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation
        for i, condition in enumerate(self.conditions):
            if condition == "control":
                pass
            """
            elif condition == 'cooperative':
                x[C.n] = 2
                x[C.KI] = 18
                x[C.K1] = 50
                x[C.K2] = 40
                x[C.K3] = 100
                x[C.K4] = 100
                x[C.K5] = 100
                x[C.K6] = 100
                x[C.K7] = 100
                x[C.K8] = 100
                x[C.K9] = 100
                x[C.K10] = 100
                x[C.V9] = 1.25
                x[C.V10] = 1.25
            """

            sol = self._solveode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[observables.index("biphosphorylated_MAPK"), :, i] = sol.y[V.MAPK_PP, :]
                self.simulations[observables.index("unphosphorylated_MAPK"), :, i] = sol.y[V.MAPK, :]

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
        Solve a system of ordinary differential equations.

        Parameters
        ----------
        diffeq : callable f(t, y, *x)
            Right-hand side of the differential equation.

        y0 : array
            Initial condition on y (can be a vector).

        t : array
            A sequence of time points for which to solve for y.

        f_params : tuple
            Model parameters.

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
        # Test data
        self.experiments[observables.index("biphosphorylated_MAPK")] = {
            "control": [
                10.0,
                298.59241599,
                298.57164292,
                295.16057672,
                211.52330623,
                80.97440738,
                84.52996829,
                296.19437698,
                254.48224797,
                122.22861484,
                36.96467096,
                296.55589221,
                284.81739793,
                166.22074898,
                47.45788287,
                215.79175866,
                293.95196966,
                210.80086084,
                80.87144427,
                84.78789794,
                296.18270875,
                254.19844509,
                121.93215568,
                37.04378877,
                296.56176029,
                284.69692824,
                165.91228275,
                47.27840604,
                216.91315533,
                293.92385296,
            ]
        }
        self.experiments[observables.index("unphosphorylated_MAPK")] = {
            "control": [
                2.80000000e02,
                1.02084698e-01,
                1.04142327e-01,
                7.54905842e-01,
                5.03960433e01,
                1.60518127e02,
                1.45421258e02,
                5.24026172e-01,
                1.98605774e01,
                1.23225761e02,
                2.02633306e02,
                4.51355981e-01,
                3.98736964e00,
                8.59280677e01,
                1.92406000e02,
                1.72867943e01,
                1.05444189e00,
                5.06418164e01,
                1.60322956e02,
                1.44990726e02,
                5.26516648e-01,
                2.00396352e01,
                1.23483781e02,
                2.02532132e02,
                4.50047126e-01,
                4.03300324e00,
                8.61819825e01,
                1.92586475e02,
                1.66298172e01,
                1.06168913e00,
            ]
        }

    @staticmethod
    def get_timepoint(obs_name):
        if obs_name in observables:
            return [60 * i for i in range(0, 150, 5)]
