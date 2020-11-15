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
        # Test data
        self.experiments[observables.index('biphosphorylated_MAPK')] = {
            'control' : [ 
                10.        , 298.59241599, 298.57164292, 295.16057672,
                211.52330623,  80.97440738,  84.52996829, 296.19437698,
                254.48224797, 122.22861484,  36.96467096, 296.55589221,
                284.81739793, 166.22074898,  47.45788287, 215.79175866,
                293.95196966, 210.80086084,  80.87144427,  84.78789794,
                296.18270875, 254.19844509, 121.93215568,  37.04378877,
                296.56176029, 284.69692824, 165.91228275,  47.27840604,
                216.91315533, 293.92385296
            ]
        }
        self.experiments[observables.index('unphosphorylated_MAPK')] = {
            'control' : [
                2.80000000e+02, 1.02084698e-01, 1.04142327e-01, 7.54905842e-01,
                5.03960433e+01, 1.60518127e+02, 1.45421258e+02, 5.24026172e-01,
                1.98605774e+01, 1.23225761e+02, 2.02633306e+02, 4.51355981e-01,
                3.98736964e+00, 8.59280677e+01, 1.92406000e+02, 1.72867943e+01,
                1.05444189e+00, 5.06418164e+01, 1.60322956e+02, 1.44990726e+02,
                5.26516648e-01, 2.00396352e+01, 1.23483781e+02, 2.02532132e+02,
                4.50047126e-01, 4.03300324e+00, 8.61819825e+01, 1.92586475e+02,
                1.66298172e+01, 1.06168913e+00
            ]
        }

    @staticmethod
    def get_timepoint(obs_name):
        if obs_name in observables:
            return [60*i for i in range(0, 150, 5)]