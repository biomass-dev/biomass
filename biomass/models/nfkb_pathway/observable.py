import numpy as np

from biomass.dynamics.solver import solve_ode

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

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[observables.index("nuclear_IkBa"), :, i] = x[C.Vnuc] * (
                    sol.y[V.nNfkIkb, :] + sol.y[V.nIkb, :]
                )
                self.simulations[observables.index("nuclear_NFkB"), :, i] = x[C.Vnuc] * (
                    sol.y[V.pnNfk, :] + sol.y[V.nNfk, :] + sol.y[V.nNfkIkb, :]
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
