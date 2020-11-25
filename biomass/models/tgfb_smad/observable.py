import numpy as np
from scipy.integrate import ode

from .name2idx import C, V
from .set_model import DifferentialEquation


observables = [
    'Ski', 'Skil', 'Dnmt3a',
    'Sox4', 'Jun', 'Smad7',
    'Klf10', 'Bmp4', 'Cxcl15',
    'Dusp5', 'Tgfa', 'Pdk4',
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
    
    t = range(600 + 1) # min

    # Experimental conditions
    conditions = [
        'WT',
        'Smad2OE',
        'Smad3OE',
        'Smad4OE'
    ]

    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation

        for gene_name in observables:
            for i, condition in enumerate(self.conditions):
                if condition == 'WT':
                    y0[V.S2] = x[C.S2tot]
                    y0[V.S3] = x[C.S3tot]
                    y0[V.S4] = x[C.S4tot]
                elif condition == 'Smad2OE':
                    y0[V.S2] = 2*x[C.S2tot]
                    y0[V.S3] = x[C.S3tot]
                    y0[V.S4] = x[C.S4tot]
                elif condition == 'Smad3OE':
                    y0[V.S2] = x[C.S2tot]
                    y0[V.S3] = 16*x[C.S3tot]
                    y0[V.S4] = x[C.S4tot]
                elif condition == 'Smad4OE':
                    y0[V.S2] = x[C.S2tot]
                    y0[V.S3] = x[C.S3tot]
                    y0[V.S4] = 3*x[C.S3tot]
                
                if gene_name == 'Ski':
                    x[C.gene_turn] = x[C.Ski_turn]
                    x[C.gene_act1] = x[C.Ski_act1]
                    x[C.gene_act2] = x[C.Ski_act2]
                    x[C.gene_act3] = x[C.Ski_act3]
                    x[C.gene_inh1] = x[C.Ski_inh1]
                    x[C.gene_inh2] = x[C.Ski_inh2]
                    x[C.gene_inh3] = x[C.Ski_inh3]
                elif gene_name == 'Skil':
                    x[C.gene_turn] = x[C.Skil_turn]
                    x[C.gene_act1] = x[C.Skil_act1]
                    x[C.gene_act2] = x[C.Skil_act2]
                    x[C.gene_act3] = x[C.Skil_act3]
                    x[C.gene_inh1] = x[C.Skil_inh1]
                    x[C.gene_inh2] = x[C.Skil_inh2]
                    x[C.gene_inh3] = x[C.Skil_inh3]
                elif gene_name == 'Dnmt3a':
                    x[C.gene_turn] = x[C.Dnmt3a_turn]
                    x[C.gene_act1] = x[C.Dnmt3a_act1]
                    x[C.gene_act2] = x[C.Dnmt3a_act2]
                    x[C.gene_act3] = x[C.Dnmt3a_act3]
                    x[C.gene_inh1] = x[C.Dnmt3a_inh1]
                    x[C.gene_inh2] = x[C.Dnmt3a_inh2]
                    x[C.gene_inh3] = x[C.Dnmt3a_inh3]
                elif gene_name == 'Sox4':
                    x[C.gene_turn] = x[C.Sox4_turn]
                    x[C.gene_act1] = x[C.Sox4_act1]
                    x[C.gene_act2] = x[C.Sox4_act2]
                    x[C.gene_act3] = x[C.Sox4_act3]
                    x[C.gene_inh1] = x[C.Sox4_inh1]
                    x[C.gene_inh2] = x[C.Sox4_inh2]
                    x[C.gene_inh3] = x[C.Sox4_inh3]
                elif gene_name == 'Jun':
                    x[C.gene_turn] = x[C.Jun_turn]
                    x[C.gene_act1] = x[C.Jun_act1]
                    x[C.gene_act2] = x[C.Jun_act2]
                    x[C.gene_act3] = x[C.Jun_act3]
                    x[C.gene_inh1] = x[C.Jun_inh1]
                    x[C.gene_inh2] = x[C.Jun_inh2]
                    x[C.gene_inh3] = x[C.Jun_inh3]
                elif gene_name == 'Smad7':
                    x[C.gene_turn] = x[C.Smad7_turn]
                    x[C.gene_act1] = x[C.Smad7_act1]
                    x[C.gene_act2] = x[C.Smad7_act2]
                    x[C.gene_act3] = x[C.Smad7_act3]
                    x[C.gene_inh1] = x[C.Smad7_inh1]
                    x[C.gene_inh2] = x[C.Smad7_inh2]
                    x[C.gene_inh3] = x[C.Smad7_inh3]
                elif gene_name == 'Klf10':
                    x[C.gene_turn] = x[C.Klf10_turn]
                    x[C.gene_act1] = x[C.Klf10_act1]
                    x[C.gene_act2] = x[C.Klf10_act2]
                    x[C.gene_act3] = x[C.Klf10_act3]
                    x[C.gene_inh1] = x[C.Klf10_inh1]
                    x[C.gene_inh2] = x[C.Klf10_inh2]
                    x[C.gene_inh3] = x[C.Klf10_inh3]
                elif gene_name == 'Bmp4':
                    x[C.gene_turn] = x[C.Bmp4_turn]
                    x[C.gene_act1] = x[C.Bmp4_act1]
                    x[C.gene_act2] = x[C.Bmp4_act2]
                    x[C.gene_act3] = x[C.Bmp4_act3]
                    x[C.gene_inh1] = x[C.Bmp4_inh1]
                    x[C.gene_inh2] = x[C.Bmp4_inh2]
                    x[C.gene_inh3] = x[C.Bmp4_inh3]
                elif gene_name == 'Cxcl15':
                    x[C.gene_turn] = x[C.Cxcl15_turn]
                    x[C.gene_act1] = x[C.Cxcl15_act1]
                    x[C.gene_act2] = x[C.Cxcl15_act2]
                    x[C.gene_act3] = x[C.Cxcl15_act3]
                    x[C.gene_inh1] = x[C.Cxcl15_inh1]
                    x[C.gene_inh2] = x[C.Cxcl15_inh2]
                    x[C.gene_inh3] = x[C.Cxcl15_inh3]
                elif gene_name == 'Dusp5':
                    x[C.gene_turn] = x[C.Dusp5_turn]
                    x[C.gene_act1] = x[C.Dusp5_act1]
                    x[C.gene_act2] = x[C.Dusp5_act2]
                    x[C.gene_act3] = x[C.Dusp5_act3]
                    x[C.gene_inh1] = x[C.Dusp5_inh1]
                    x[C.gene_inh2] = x[C.Dusp5_inh2]
                    x[C.gene_inh3] = x[C.Dusp5_inh3]
                elif gene_name == 'Tgfa':
                    x[C.gene_turn] = x[C.Tgfa_turn]
                    x[C.gene_act1] = x[C.Tgfa_act1]
                    x[C.gene_act2] = x[C.Tgfa_act2]
                    x[C.gene_act3] = x[C.Tgfa_act3]
                    x[C.gene_inh1] = x[C.Tgfa_inh1]
                    x[C.gene_inh2] = x[C.Tgfa_inh2]
                    x[C.gene_inh3] = x[C.Tgfa_inh3]
                elif gene_name == 'Pdk4':
                    x[C.gene_turn] = x[C.Pdk4_turn]
                    x[C.gene_act1] = x[C.Pdk4_act1]
                    x[C.gene_act2] = x[C.Pdk4_act2]
                    x[C.gene_act3] = x[C.Pdk4_act3]
                    x[C.gene_inh1] = x[C.Pdk4_inh1]
                    x[C.gene_inh2] = x[C.Pdk4_inh2]
                    x[C.gene_inh3] = x[C.Pdk4_inh3]
                
                (T, Y) = self._solveode(self.diffeq, y0, self.t, tuple(x))

                if T[-1] < self.t[-1]:
                    return False
                else:
                    self.simulations[observables.index(gene_name), :, i] = (
                        np.log2(Y[:, V.gene])
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