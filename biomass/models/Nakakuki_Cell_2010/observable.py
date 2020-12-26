import numpy as np
from scipy.integrate import solve_ivp
from typing import List, Callable, Optional

from .name2idx import C, V
from .set_model import DifferentialEquation


observables = [
    "Phosphorylated_MEKc",
    "Phosphorylated_ERKc",
    "Phosphorylated_RSKw",
    "Phosphorylated_CREBw",
    "dusp_mRNA",
    "cfos_mRNA",
    "cFos_Protein",
    "Phosphorylated_cFos",
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
        for observable in observables:
            self.normalization[observable] = {"timepoint": None, "condition": []}

    t = range(5401)  # 0, 1, 2, ..., 5400 (Unit: sec.)

    # Experimental conditions
    conditions = ["EGF", "HRG"]

    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation
        # get steady state
        x[C.Ligand] = x[C.no_ligand]  # No ligand
        y0 = self._get_steady_state(self.diffeq, y0, tuple(x))
        if not y0:
            return False
        # add ligand
        for i, condition in enumerate(self.conditions):
            if condition == "EGF":
                x[C.Ligand] = x[C.EGF]
            elif condition == "HRG":
                x[C.Ligand] = x[C.HRG]

            sol = self._solveode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[observables.index("Phosphorylated_MEKc"), :, i] = sol.y[V.ppMEKc, :]
                self.simulations[observables.index("Phosphorylated_ERKc"), :, i] = (
                    sol.y[V.pERKc, :] + sol.y[V.ppERKc, :]
                )
                self.simulations[observables.index("Phosphorylated_RSKw"), :, i] = sol.y[V.pRSKc, :] + sol.y[
                    V.pRSKn, :
                ] * (x[C.Vn] / x[C.Vc])
                self.simulations[observables.index("Phosphorylated_CREBw"), :, i] = sol.y[V.pCREBn, :] * (
                    x[C.Vn] / x[C.Vc]
                )
                self.simulations[observables.index("dusp_mRNA"), :, i] = sol.y[V.duspmRNAc, :]
                self.simulations[observables.index("cfos_mRNA"), :, i] = sol.y[V.cfosmRNAc, :]
                self.simulations[observables.index("cFos_Protein"), :, i] = (
                    (sol.y[V.pcFOSn, :] + sol.y[V.cFOSn, :]) * (x[C.Vn] / x[C.Vc])
                    + sol.y[V.cFOSc, :]
                    + sol.y[V.pcFOSc, :]
                )
                self.simulations[observables.index("Phosphorylated_cFos"), :, i] = (
                    sol.y[V.pcFOSn, :] * (x[C.Vn] / x[C.Vc]) + sol.y[V.pcFOSc, :]
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
        self.experiments[observables.index("Phosphorylated_MEKc")] = {
            "EGF": [0.000, 0.773, 0.439, 0.252, 0.130, 0.087, 0.080, 0.066],
            "HRG": [0.000, 0.865, 1.000, 0.837, 0.884, 0.920, 0.875, 0.789],
        }
        self.error_bars[observables.index("Phosphorylated_MEKc")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0.000, 0.030, 0.048, 0.009, 0.009, 0.017, 0.012, 0.008]],
            "HRG": [sd / np.sqrt(3) for sd in [0.000, 0.041, 0.000, 0.051, 0.058, 0.097, 0.157, 0.136]],
        }

        self.experiments[observables.index("Phosphorylated_ERKc")] = {
            "EGF": [0.000, 0.867, 0.799, 0.494, 0.313, 0.266, 0.200, 0.194],
            "HRG": [0.000, 0.848, 1.000, 0.971, 0.950, 0.812, 0.747, 0.595],
        }
        self.error_bars[observables.index("Phosphorylated_ERKc")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0.000, 0.137, 0.188, 0.126, 0.096, 0.087, 0.056, 0.012]],
            "HRG": [sd / np.sqrt(3) for sd in [0.000, 0.120, 0.000, 0.037, 0.088, 0.019, 0.093, 0.075]],
        }

        self.experiments[observables.index("Phosphorylated_RSKw")] = {
            "EGF": [0, 0.814, 0.812, 0.450, 0.151, 0.059, 0.038, 0.030],
            "HRG": [0, 0.953, 1.000, 0.844, 0.935, 0.868, 0.779, 0.558],
        }
        self.error_bars[observables.index("Phosphorylated_RSKw")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0, 0.064, 0.194, 0.030, 0.027, 0.031, 0.043, 0.051]],
            "HRG": [sd / np.sqrt(3) for sd in [0, 0.230, 0.118, 0.058, 0.041, 0.076, 0.090, 0.077]],
        }

        self.experiments[observables.index("Phosphorylated_cFos")] = {
            "EGF": [0, 0.060, 0.109, 0.083, 0.068, 0.049, 0.027, 0.017],
            "HRG": [0, 0.145, 0.177, 0.158, 0.598, 1.000, 0.852, 0.431],
        }
        self.error_bars[observables.index("Phosphorylated_cFos")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0, 0.003, 0.021, 0.013, 0.016, 0.007, 0.003, 0.002]],
            "HRG": [sd / np.sqrt(3) for sd in [0, 0.010, 0.013, 0.001, 0.014, 0.000, 0.077, 0.047]],
        }

        # ----------------------------------------------------------------------

        self.experiments[observables.index("Phosphorylated_CREBw")] = {
            "EGF": [0, 0.446, 0.030, 0.000, 0.000],
            "HRG": [0, 1.000, 0.668, 0.460, 0.340],
        }
        self.error_bars[observables.index("Phosphorylated_CREBw")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]],
            "HRG": [sd / np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]],
        }
        # ----------------------------------------------------------------------

        self.experiments[observables.index("cfos_mRNA")] = {
            "EGF": [0, 0.181, 0.476, 0.518, 0.174, 0.026, 0.000],
            "HRG": [0, 0.353, 0.861, 1.000, 0.637, 0.300, 0.059],
        }
        self.error_bars[observables.index("cfos_mRNA")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0.017, 0.004, 0.044, 0.004, 0.023, 0.007, 0.008]],
            "HRG": [sd / np.sqrt(3) for sd in [0.017, 0.006, 0.065, 0.044, 0.087, 0.023, 0.001]],
        }
        # ----------------------------------------------------------------------

        self.experiments[observables.index("cFos_Protein")] = {
            "EGF": [0, 0.078, 0.216, 0.240, 0.320, 0.235],
            "HRG": [0, 0.089, 0.552, 0.861, 1.000, 0.698],
        }
        self.error_bars[observables.index("cFos_Protein")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0, 0.036, 0.028, 0.056, 0.071, 0.048]],
            "HRG": [sd / np.sqrt(3) for sd in [0, 0.021, 0.042, 0.063, 0.000, 0.047]],
        }

        self.experiments[observables.index("dusp_mRNA")] = {
            "EGF": [0.000, 0.177, 0.331, 0.214, 0.177, 0.231],
            "HRG": [0.000, 0.221, 0.750, 1.000, 0.960, 0.934],
        }
        self.error_bars[observables.index("dusp_mRNA")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0.033, 0.060, 0.061, 0.032, 0.068, 0.050]],
            "HRG": [sd / np.sqrt(3) for sd in [0.027, 0.059, 0.094, 0.124, 0.113, 0.108]],
        }

    @staticmethod
    def get_timepoint(obs_name):
        if obs_name in ["Phosphorylated_MEKc", "Phosphorylated_ERKc", "Phosphorylated_RSKw", "Phosphorylated_cFos"]:
            return [0, 300, 600, 900, 1800, 2700, 3600, 5400]  # (Unit: sec.)

        elif obs_name == "Phosphorylated_CREBw":
            return [0, 600, 1800, 3600, 5400]

        elif obs_name == "cfos_mRNA":
            return [0, 600, 1200, 1800, 2700, 3600, 5400]

        elif obs_name in ["cFos_Protein", "dusp_mRNA"]:
            return [0, 900, 1800, 2700, 3600, 5400]