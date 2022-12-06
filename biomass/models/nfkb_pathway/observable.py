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
        Names of model observables.

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
        self.obs_names: list = [
            "nuclear_IkBa",
            "nuclear_NFkB",
        ]
        self.t: range = range(200 + 1)
        self.conditions: list = ["TNFa", "TNFa_DCF"]
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.conditions), len(self.t))
        )
        self.normalization = {}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

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

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index("nuclear_IkBa"), i] = x[C.Vnuc] * (
                    sol.y[V.nNfkIkb] + sol.y[V.nIkb]
                )
                self.simulations[self.obs_names.index("nuclear_NFkB"), i] = x[C.Vnuc] * (
                    sol.y[V.pnNfk] + sol.y[V.nNfk] + sol.y[V.nNfkIkb]
                )
        return None

    def set_data(self):
        pass

    def get_timepoint(self, obs_name):
        if obs_name in self.obs_names:
            return []
        assert False
