import numpy as np

from biomass.dynamics.solver import solve_ode

from .name2idx import C, V
from .ode import DifferentialEquation

class Observable(DifferentialEquation):

    def __init__(self):
        super(Observable, self).__init__(pertubation={})
        self.obs_names: list = [
            "Complex",
            "Substrate",
            "Product",
        ]
        self.t: range = range(100)
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.t))
        )
        self.normalization: dict = {}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

    def simulate(self, x, y0, _pertubation={}):
        if _pertubation:
            self.pertubation = _pertubation
        sol = solve_ode(self.diffeq, y0, self.t, tuple(x))
        if sol is None:
            return False
        else:
            self.simulations[self.obs_names.index("Complex")] = sol.y[V.ES]
            self.simulations[self.obs_names.index("Substrate")] = sol.y[V.S]
            self.simulations[self.obs_names.index("Product")] = sol.y[V.P]

    def set_data(self):
        pass

    def get_timepoint(self, obs_name):
        if obs_name in self.obs_names:
            return []
