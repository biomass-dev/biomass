from typing import List

import numpy as np

from biomass.dynamics.solver import get_steady_state, solve_ode

from .name2idx import C, V
from .ode import DifferentialEquation


class Observable(DifferentialEquation):
    """
    Correlating model simulations and experimental measurements.

    Attributes
    ----------
    obs_names : list of strings
        Names of self.obs_names.

    t : range
        Simulation time span.

    conditions : list of strings
        Experimental conditions.

    simulations : numpy.ndarray
        The numpy array to store simulation results.

    normalization : nested dict
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If :obj:`None`, the maximum value will be used for normalization.

        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in ``sim.conditions`` will be used.

    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """

    def __init__(self):
        super(Observable, self).__init__(perturbation={})
        self.obs_names = [
            "Phosphorylated_MEKc",
            "Phosphorylated_ERKc",
            "Phosphorylated_RSKw",
            "Phosphorylated_CREBw",
            "dusp_mRNA",
            "cfos_mRNA",
            "cFos_Protein",
            "Phosphorylated_cFos",
        ]
        self.t: range = range(5401)
        self.conditions: list = ["EGF", "HRG"]
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.conditions), len(self.t))
        )
        self.normalization: dict = {}
        for observable in self.obs_names:
            self.normalization[observable] = {"timepoint": None, "condition": []}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

    def simulate(self, x, y0, _perturbation=None):
        if _perturbation is not None:
            self.perturbation = _perturbation
        # get steady state
        x[C.Ligand] = x[C.no_ligand]  # No ligand
        y0 = get_steady_state(self.diffeq, y0, tuple(x), integrator="zvode")
        if not y0:
            return False
        # add ligand
        for i, condition in enumerate(self.conditions):
            if condition == "EGF":
                x[C.Ligand] = x[C.EGF]
            elif condition == "HRG":
                x[C.Ligand] = x[C.HRG]

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x), method="BDF")

            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index("Phosphorylated_MEKc"), i] = sol.y[V.ppMEKc]
                self.simulations[self.obs_names.index("Phosphorylated_ERKc"), i] = (
                    sol.y[V.pERKc] + sol.y[V.ppERKc]
                )
                self.simulations[self.obs_names.index("Phosphorylated_RSKw"), i] = sol.y[
                    V.pRSKc
                ] + sol.y[V.pRSKn] * (x[C.Vn] / x[C.Vc])
                self.simulations[self.obs_names.index("Phosphorylated_CREBw"), i] = sol.y[
                    V.pCREBn
                ] * (x[C.Vn] / x[C.Vc])
                self.simulations[self.obs_names.index("dusp_mRNA"), i] = sol.y[V.duspmRNAc]
                self.simulations[self.obs_names.index("cfos_mRNA"), i] = sol.y[V.cfosmRNAc]
                self.simulations[self.obs_names.index("cFos_Protein"), i] = (
                    (sol.y[V.pcFOSn] + sol.y[V.cFOSn]) * (x[C.Vn] / x[C.Vc])
                    + sol.y[V.cFOSc]
                    + sol.y[V.pcFOSc]
                )
                self.simulations[self.obs_names.index("Phosphorylated_cFos"), i] = (
                    sol.y[V.pcFOSn] * (x[C.Vn] / x[C.Vc]) + sol.y[V.pcFOSc]
                )
        return None

    def set_data(self):
        self.experiments[self.obs_names.index("Phosphorylated_MEKc")] = {
            "EGF": [0.000, 0.773, 0.439, 0.252, 0.130, 0.087, 0.080, 0.066],
            "HRG": [0.000, 0.865, 1.000, 0.837, 0.884, 0.920, 0.875, 0.789],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_MEKc")] = {
            "EGF": [
                sd / np.sqrt(3) for sd in [0.000, 0.030, 0.048, 0.009, 0.009, 0.017, 0.012, 0.008]
            ],
            "HRG": [
                sd / np.sqrt(3) for sd in [0.000, 0.041, 0.000, 0.051, 0.058, 0.097, 0.157, 0.136]
            ],
        }

        self.experiments[self.obs_names.index("Phosphorylated_ERKc")] = {
            "EGF": [0.000, 0.867, 0.799, 0.494, 0.313, 0.266, 0.200, 0.194],
            "HRG": [0.000, 0.848, 1.000, 0.971, 0.950, 0.812, 0.747, 0.595],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_ERKc")] = {
            "EGF": [
                sd / np.sqrt(3) for sd in [0.000, 0.137, 0.188, 0.126, 0.096, 0.087, 0.056, 0.012]
            ],
            "HRG": [
                sd / np.sqrt(3) for sd in [0.000, 0.120, 0.000, 0.037, 0.088, 0.019, 0.093, 0.075]
            ],
        }

        self.experiments[self.obs_names.index("Phosphorylated_RSKw")] = {
            "EGF": [0, 0.814, 0.812, 0.450, 0.151, 0.059, 0.038, 0.030],
            "HRG": [0, 0.953, 1.000, 0.844, 0.935, 0.868, 0.779, 0.558],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_RSKw")] = {
            "EGF": [
                sd / np.sqrt(3) for sd in [0, 0.064, 0.194, 0.030, 0.027, 0.031, 0.043, 0.051]
            ],
            "HRG": [
                sd / np.sqrt(3) for sd in [0, 0.230, 0.118, 0.058, 0.041, 0.076, 0.090, 0.077]
            ],
        }

        self.experiments[self.obs_names.index("Phosphorylated_cFos")] = {
            "EGF": [0, 0.060, 0.109, 0.083, 0.068, 0.049, 0.027, 0.017],
            "HRG": [0, 0.145, 0.177, 0.158, 0.598, 1.000, 0.852, 0.431],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_cFos")] = {
            "EGF": [
                sd / np.sqrt(3) for sd in [0, 0.003, 0.021, 0.013, 0.016, 0.007, 0.003, 0.002]
            ],
            "HRG": [
                sd / np.sqrt(3) for sd in [0, 0.010, 0.013, 0.001, 0.014, 0.000, 0.077, 0.047]
            ],
        }

        # ----------------------------------------------------------------------

        self.experiments[self.obs_names.index("Phosphorylated_CREBw")] = {
            "EGF": [0, 0.446, 0.030, 0.000, 0.000],
            "HRG": [0, 1.000, 0.668, 0.460, 0.340],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_CREBw")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]],
            "HRG": [sd / np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]],
        }
        # ----------------------------------------------------------------------

        self.experiments[self.obs_names.index("cfos_mRNA")] = {
            "EGF": [0, 0.181, 0.476, 0.518, 0.174, 0.026, 0.000],
            "HRG": [0, 0.353, 0.861, 1.000, 0.637, 0.300, 0.059],
        }
        self.error_bars[self.obs_names.index("cfos_mRNA")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0.017, 0.004, 0.044, 0.004, 0.023, 0.007, 0.008]],
            "HRG": [sd / np.sqrt(3) for sd in [0.017, 0.006, 0.065, 0.044, 0.087, 0.023, 0.001]],
        }
        # ----------------------------------------------------------------------

        self.experiments[self.obs_names.index("cFos_Protein")] = {
            "EGF": [0, 0.078, 0.216, 0.240, 0.320, 0.235],
            "HRG": [0, 0.089, 0.552, 0.861, 1.000, 0.698],
        }
        self.error_bars[self.obs_names.index("cFos_Protein")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0, 0.036, 0.028, 0.056, 0.071, 0.048]],
            "HRG": [sd / np.sqrt(3) for sd in [0, 0.021, 0.042, 0.063, 0.000, 0.047]],
        }

        self.experiments[self.obs_names.index("dusp_mRNA")] = {
            "EGF": [0.000, 0.177, 0.331, 0.214, 0.177, 0.231],
            "HRG": [0.000, 0.221, 0.750, 1.000, 0.960, 0.934],
        }
        self.error_bars[self.obs_names.index("dusp_mRNA")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0.033, 0.060, 0.061, 0.032, 0.068, 0.050]],
            "HRG": [sd / np.sqrt(3) for sd in [0.027, 0.059, 0.094, 0.124, 0.113, 0.108]],
        }

    @staticmethod
    def get_timepoint(obs_name) -> List[int]:
        if obs_name in [
            "Phosphorylated_MEKc",
            "Phosphorylated_ERKc",
            "Phosphorylated_RSKw",
            "Phosphorylated_cFos",
        ]:
            return [0, 300, 600, 900, 1800, 2700, 3600, 5400]  # (Unit: sec.)
        elif obs_name == "Phosphorylated_CREBw":
            return [0, 600, 1800, 3600, 5400]
        elif obs_name == "cfos_mRNA":
            return [0, 600, 1200, 1800, 2700, 3600, 5400]
        elif obs_name in ["cFos_Protein", "dusp_mRNA"]:
            return [0, 900, 1800, 2700, 3600, 5400]
        assert False
