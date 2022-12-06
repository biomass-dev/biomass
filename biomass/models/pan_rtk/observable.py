import numpy as np

from biomass.dynamics.solver import solve_ode

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
            "pEGFR_au",
            "pErbB2_au",
            "pErbB3_au",
            "pERK_au",
            "pAKT_au",
            "pS6_au",
        ]
        self.t: range = range(241)
        self.conditions: list = ["EGF0nM", "EGF0_156nM", "EGF0_625nM", "EGF2_5nM", "EGF10nM"]
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.conditions), len(self.t))
        )
        self.normalization: dict = {}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

    def simulate(self, x, y0, _perturbation=None):
        if _perturbation is not None:
            self.perturbation = _perturbation
        # add ligand
        for i, condition in enumerate(self.conditions):
            if condition == "EGF0nM":
                y0[V.dose_EGF] = 0
                y0[V.dose_HGF] = 0
                y0[V.dose_IGF1] = 0
                y0[V.dose_HRG] = 0
            elif condition == "EGF0_156nM":
                y0[V.dose_EGF] = 0.156 * x[C.scale_Ligand]
                y0[V.dose_IGF1] = 0
                y0[V.dose_HRG] = 0
            elif condition == "EGF0_625nM":
                y0[V.dose_EGF] = 0.625 * x[C.scale_Ligand]
                y0[V.dose_HGF] = 0
                y0[V.dose_IGF1] = 0
                y0[V.dose_HRG] = 0
            elif condition == "EGF2_5nM":
                y0[V.dose_EGF] = 2.5 * x[C.scale_Ligand]
                y0[V.dose_HGF] = 0
                y0[V.dose_IGF1] = 0
                y0[V.dose_HRG] = 0
            elif condition == "EGF10nM":
                y0[V.dose_EGF] = 10 * x[C.scale_Ligand]
                y0[V.dose_HGF] = 0
                y0[V.dose_IGF1] = 0
                y0[V.dose_HRG] = 0

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index("pEGFR_au"), i] = np.log10(
                    x[C.offset_pEGFR_CelllineH322M]
                    + x[C.scale_pEGFR_CelllineH322M]
                    * (
                        2 * sol.y[V.pEGFRd]
                        + 2 * sol.y[V.pEGFRi]
                        + 2 * sol.y[V.pEGFRi_ph]
                        + sol.y[V.pErbB12]
                        + sol.y[V.pErbB12i]
                        + sol.y[V.pErbB12i_ph]
                        + sol.y[V.pErbB13]
                        + sol.y[V.pErbB13i]
                        + sol.y[V.pErbB13i_ph]
                        + sol.y[V.pMetEGFR]
                        + sol.y[V.pMetEGFRi]
                        + sol.y[V.pMetEGFRi_ph]
                    )
                )
                self.simulations[self.obs_names.index("pErbB2_au"), i] = np.log10(
                    x[C.offset_pErbB2_CelllineH322M]
                    + x[C.scale_pErbB2_CelllineH322M]
                    * (
                        sol.y[V.pErbB12]
                        + sol.y[V.pErbB12i]
                        + sol.y[V.pErbB12i_ph]
                        + 2 * sol.y[V.pErbB2]
                        + 2 * sol.y[V.pErbB2i]
                        + 2 * sol.y[V.pErbB2i_ph]
                        + sol.y[V.pErbB32]
                        + sol.y[V.pErbB32i]
                        + sol.y[V.pErbB32i_ph]
                    )
                )
                self.simulations[self.obs_names.index("pErbB3_au"), i] = np.log10(
                    x[C.offset_pErbB3_CelllineH322M]
                    + x[C.scale_pErbB3_CelllineH322M]
                    * (
                        sol.y[V.pErbB13]
                        + sol.y[V.pErbB13i]
                        + sol.y[V.pErbB13i_ph]
                        + sol.y[V.pErbB32]
                        + sol.y[V.pErbB32i]
                        + sol.y[V.pErbB32i_ph]
                        + 2 * sol.y[V.pErbB3d]
                        + 2 * sol.y[V.pErbB3i]
                        + 2 * sol.y[V.pErbB3i_ph]
                        + sol.y[V.pMetErbB3]
                        + sol.y[V.pMetErbB3i]
                        + sol.y[V.pMetErbB3i_ph]
                    )
                )
                self.simulations[self.obs_names.index("pERK_au"), i] = np.log10(
                    x[C.offset_pERK_CelllineH322M]
                    + x[C.scale_pERK_CelllineH322M] * (sol.y[V.pERK])
                )
                self.simulations[self.obs_names.index("pAKT_au"), i] = np.log10(
                    x[C.offset_pAKT_CelllineH322M] + x[C.scale_pAKT_CelllineH322M] * sol.y[V.pAKT]
                )
                self.simulations[self.obs_names.index("pS6_au"), i] = np.log10(
                    x[C.offset_pS6_CelllineH322M] + x[C.scale_pS6_CelllineH322M] * sol.y[V.pS6]
                )
        return None

    def set_data(self):
        pass

    def get_timepoint(self, obs_name):
        if obs_name in self.obs_names:
            return []
        assert False
