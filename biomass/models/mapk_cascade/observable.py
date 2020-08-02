import numpy as np
from scipy.integrate import ode

from .name2idx import C, V
from .set_model import DifferentialEquation


observables = [
    'biphosphorylated_MAPK',
    'unphosphorylated_MAPK',
]

class NumericalSimulation(DifferentialEquation):
    """ Simulate a model using scipy.integrate.ode

    Attributes
    ----------
    normalization : bool
        if True, simulation results in each observable are divided by their 
        maximum values.

    """
    def __init__(self):
        super().__init__(perturbation={})
        self.normalization = False

    t = range(150*60+1)

    # Experimental conditions
    conditions = ['control']

    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation
        for i, condition in enumerate(self.conditions):
            if condition == 'control':
                pass
            '''
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
            '''
            
            (T, Y) = self._solveode(self.diffeq, y0, self.t, tuple(x))

            if T[-1] < self.t[-1]:
                return False
            else:
                self.simulations[observables.index('biphosphorylated_MAPK'), :, i] = (
                    Y[:, V.MAPK_PP]
                )
                self.simulations[observables.index('unphosphorylated_MAPK'), :, i] = (
                    Y[:, V.MAPK]
                )

    def _solveode(self, diffeq, y0, tspan, args):
        dt = (self.t[-1] - self.t[0]) / (len(self.t) - 1)
        sol = ode(diffeq)
        sol.set_integrator(
            'vode', method='bdf', with_jacobian=True,
            atol=1e-9, rtol=1e-9, min_step=1e-8
        )
        sol.set_initial_value(y0, tspan[0])
        sol.set_f_params(args)

        T = [tspan[0]]
        Y = [y0]

        while sol.successful() and sol.t < tspan[-1]:
            sol.integrate(sol.t+dt)
            T.append(sol.t)
            Y.append(sol.y)

        return np.array(T), np.array(Y)
    
    def _get_steady_state(self, diffeq, y0, args, eps=1e-6):
        """
        Run until a time t for which the maximal absolutevalue of the 
        regularized relative derivative was smaller than eps.
        """
        while True:
            (T, Y) = self._solveode(diffeq, y0, range(2), args)
            if T[-1] < 1:
                return None
            elif np.max(np.abs((Y[-1, :] - y0) / (np.array(y0) + eps))) < eps:
                break
            else:
                y0 = Y[-1, :].tolist()

        return y0

class ExperimentalData(object):
    def __init__(self):
        pass

    experiments = [None] * len(observables)
    error_bar = [None] * len(observables)

    def get_timepoint(self, obs_idx):
        '''
        return list(map(int, exp_t))
        '''
        pass