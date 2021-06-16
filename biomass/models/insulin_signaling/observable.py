import numpy as np

from biomass.dynamics.solver import solve_ode

from .name2idx import C, V
from .set_model import DifferentialEquation


class Observable(DifferentialEquation):
    """
    Correlating model simulations and experimental measurements.

    Attributes
    ----------
    obs_names : list of strings
        Names of model observables.

    t : range
        Simulation time span.

    conditions : list of strings
        Expetimental conditions.

    simulations : numpy.ndarray
        The numpy array to store simulation results.

    normalization : nested dict
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If None, the maximum value will be used for normalization.

        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in sim.conditions will be used.

    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """

    def __init__(self):
        super(Observable, self).__init__(perturbation={})
        self.obs_names: list = [
            "pAKT",
            "pS6K",
            "pGSK3B",
            "G6Pase",
        ]
        self.t: range = range(480 + 1)  # min
        self.conditions: list = [
            "Insulin 0.01",
            "Insulin 0.03",
            "Insulin 0.1",
            "Insulin 0.3",
            "Insulin 1.0",
        ]
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.t), len(self.conditions))
        )
        self.normalization: dict = {}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation

        y0[V.Ins] = 0.01  # 0.01 nM of insulin during starvation
        sol = solve_ode(self.diffeq, y0, range(2400 + 1), tuple(x))
        if sol is None:
            return False
        else:
            y0 = sol.y[:, -1].tolist()

        for i, condition in enumerate(self.conditions):
            if condition == "Insulin 0.01":
                y0[V.Ins] = 0.01
            elif condition == "Insulin 0.03":
                y0[V.Ins] = 0.03
            elif condition == "Insulin 0.1":
                y0[V.Ins] = 0.1
            elif condition == "Insulin 0.3":
                y0[V.Ins] = 0.3
            elif condition == "Insulin 1.0":
                y0[V.Ins] = 1.0

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index("pAKT"), :, i] = sol.y[V.pAKT, :]
                self.simulations[self.obs_names.index("pS6K"), :, i] = (
                    sol.y[V.pS6K, :] * 83.8672192461257
                )
                self.simulations[self.obs_names.index("pGSK3B"), :, i] = (
                    sol.y[V.pGSK3B, :] * 0.111097316860158
                )
                self.simulations[self.obs_names.index("G6Pase"), :, i] = (
                    sol.y[V.G6Pase, :] * 0.0363622452066626
                )

    def set_data(self):
        pass

    def get_timepoint(self, obs_name):
        if obs_name in self.obs_names:
            return []
