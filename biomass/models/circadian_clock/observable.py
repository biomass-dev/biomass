import numpy as np
from scipy.integrate import ode

from .name2idx import C, V
from .set_model import DifferentialEquation


observables = [
    'Per_mRNA',
    'Cry_mRNA',
    'Bmal1_mRNA',
]

class NumericalSimulation(DifferentialEquation):
    """ Simulate a model using scipy.integrate.ode

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

    t = range(72+1)

    # Experimental conditions
    conditions = ['DD']

    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation
        for i, condition in enumerate(self.conditions):
            if condition == 'DD':
                pass
            
            (T, Y) = self._solveode(self.diffeq, y0, self.t, tuple(x))

            if T[-1] < self.t[-1]:
                return False
            else:
                self.simulations[observables.index('Per_mRNA'), :, i] = (
                    Y[:, V.MP]
                )
                self.simulations[observables.index('Cry_mRNA'), :, i] = (
                    Y[:, V.MC]
                )
                self.simulations[observables.index('Bmal1_mRNA'), :, i] = (
                    Y[:, V.MB]
                )

    def _solveode(self, diffeq, y0, tspan, args):
        """
        Solve a system of ordinary differential equations.

        Parameters
        ----------
        diffeq : callable f(y, t, f_args)
            Right-hand side of the differential equation.

        y0 : array
            Initial condition on y (can be a vector).
        
        tspan : array
            A sequence of time points for which to solve for y.
        
        args : tuple
            Model parameters.
        
        Returns
        -------
        T, Y : tuple
            T : array, shape (len(t))
                Evaluation points.

            Y : array, shape (len(t), len(y0))
                Array containing the value of y for each desired time in t, 
                with the initial value y0 in the first row.

        """
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
        Find the steady state for the untreated condition.

        Parameters
        ----------
        diffeq : callable f(y, t, f_args)
            Right-hand side of the differential equation.

        y0 : array
            Initial condition on y (can be a vector).
        
        args : tuple
            Model parameters.
        
        eps : float (default: 1e-6)
            Run until a time t for which the maximal absolutevalue of the 
            regularized relative derivative was smaller than eps.
        
        Returns
        -------
        y0 : array
            Steady state concentrations of all species.
            
        """
        while True:
            (T, Y) = self._solveode(diffeq, y0, range(2), args)
            if T[-1] < 1 or \
                    np.max(
                        np.abs((Y[-1, :] - y0) / (np.array(y0) + eps))
                    ) < eps:
                break
            else:
                y0 = Y[-1, :].tolist()

        return [] if T[-1] < 1 else Y[-1, :].tolist()

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