import numpy as np

from biomass.dynamics.solver import solve_ode

from .name2idx import C, V
from .set_model import DifferentialEquation

observables = [
    "pAKT",
    "pS6K",
    "pGSK3B",
    "G6Pase",
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

    t = range(480 + 1)  # min

    # Experimental conditions
    conditions = ["Insulin 0.01", "Insulin 0.03", "Insulin 0.1", "Insulin 0.3", "Insulin 1.0"]

    simulations = np.empty((len(observables), len(t), len(conditions)))

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
                self.simulations[observables.index("pAKT"), :, i] = sol.y[V.pAKT, :]
                self.simulations[observables.index("pS6K"), :, i] = (
                    sol.y[V.pS6K, :] * 83.8672192461257
                )
                self.simulations[observables.index("pGSK3B"), :, i] = (
                    sol.y[V.pGSK3B, :] * 0.111097316860158
                )
                self.simulations[observables.index("G6Pase"), :, i] = (
                    sol.y[V.G6Pase, :] * 0.0363622452066626
                )


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
