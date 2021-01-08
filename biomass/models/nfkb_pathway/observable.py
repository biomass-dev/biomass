import numpy as np
from scipy.integrate import solve_ivp
from typing import List, Callable, Optional

from .name2idx import C, V
from .set_model import DifferentialEquation


observables = [
    "nuclear_IkBa",
    "nuclear_NFkB",
]


class NumericalSimulation(DifferentialEquation):
    """Simulate a model using scipy.integrate.solve_ivp

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

    t = range(200 + 1)

    # Experimental conditions
    conditions = ["TNFa", "TNFa_DCF"]

    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation
        for i, condition in enumerate(self.conditions):
            if condition == "TNFa":
                pass
            elif condition == "TNFa_DCF":
                x[C.uptake] = 1.0000
                x[C.TNF] = 1.0000
                x[C.trigger_iIkk] = 0.0195
                x[C.deact_TNFR] = 0.0010
                x[C.deact_ppIkk] = 0.1660
                x[C.deact_pnNfk] = 1000.0000
                x[C.act_Ikk_by_TNF] = 0.0347
                x[C.act_pIkk] = 0.1603
                x[C.act_Ikb_by_Ikk] = 0.1562
                x[C.act_Nfk_by_Ikk] = 0.6438
                x[C.act_Nfk_by_Ikk_complex] = 0.2816
                x[C.act_Ikb_complex] = 1.3897
                x[C.form_complex] = 2.8390
                x[C.form_complex_nuc] = 1000.0000
                x[C.ext_nNfkIkb] = 1000.0000
                x[C.Vnuc] = 1.0000
                x[C.split_NfkpIkb] = 0.0811
                x[C.split_NfkIkb] = 1.0000
                x[C.int_Nfk] = 0.0100
                x[C.int_Ikb] = 0.1226
                x[C.eta_int_pNfk] = 17.9585
                x[C.degrad_Ikb] = 0.6308
                x[C.degrad_mIkb] = 0.0053
                x[C.degrad_RnaA20] = 0.0089
                x[C.degrad_A20] = 0.0116
                x[C.prod_Ikb] = 1.0000
                x[C.prod_mIkb_by_nNfk] = 0.0020
                x[C.build_RnaA20] = 1.0000
                x[C.build_A20] = 0.0006
                x[C.shuttle_RnaA20] = 0.0119

            sol = self._solveode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[observables.index("nuclear_IkBa"), :, i] = x[C.Vnuc] * (
                    sol.y[V.nNfkIkb, :] + sol.y[V.nIkb, :]
                )
                self.simulations[observables.index("nuclear_NFkB"), :, i] = x[C.Vnuc] * (
                    sol.y[V.pnNfk, :] + sol.y[V.nNfk, :] + sol.y[V.nNfkIkb, :]
                )

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
        pass

    @staticmethod
    def get_timepoint(obs_name):
        if obs_name in observables:
            return []