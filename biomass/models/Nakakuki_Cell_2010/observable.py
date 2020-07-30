import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import ode

from .name2idx import C, V
from .set_model import DifferentialEquation


observables = [
    'Phosphorylated_MEKc',
    'Phosphorylated_ERKc',
    'Phosphorylated_RSKw',
    'Phosphorylated_CREBw',
    'dusp_mRNA',
    'cfos_mRNA',
    'cFos_Protein',
    'Phosphorylated_cFos',
]

class NumericalSimulation(DifferentialEquation):
    def __init__(self):
        super().__init__(perturbation={})
        self.normalization = True
        '''
        if True, simulation results in each observable 
        are divided by their maximum values
        '''

    t = range(5401)  # 0, 1, 2, ..., 5400 (Unit: sec.)

    # Experimental conditions
    conditions = ['EGF', 'HRG']

    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation
        # get steady state
        x[C.Ligand] = x[C.no_ligand]  # No ligand
        Y_steady_state = self._get_steady_state(self.diffeq, y0, tuple(x))
        if Y_steady_state is None:
            return False
        else:
            y0 = Y_steady_state[:]
        # add ligand
        for i, condition in enumerate(self.conditions):
            if condition == 'EGF':
                x[C.Ligand] = x[C.EGF]
            elif condition == 'HRG':
                x[C.Ligand] = x[C.HRG]

            (T, Y) = self._solveode(self.diffeq, y0, self.t, tuple(x))

            if T[-1] < self.t[-1]:
                return False
            else:
                self.simulations[observables.index('Phosphorylated_MEKc'), :, i] = (
                    Y[:, V.ppMEKc]
                )
                self.simulations[observables.index('Phosphorylated_ERKc'), :, i] = (
                    Y[:, V.pERKc] + Y[:, V.ppERKc]
                )
                self.simulations[observables.index('Phosphorylated_RSKw'), :, i] = (
                    Y[:, V.pRSKc] + Y[:, V.pRSKn] * (x[C.Vn]/x[C.Vc])
                )
                self.simulations[observables.index('Phosphorylated_CREBw'), :, i] = (
                    Y[:, V.pCREBn]*(x[C.Vn]/x[C.Vc])
                )
                self.simulations[observables.index('dusp_mRNA'), :, i] = (
                    Y[:, V.duspmRNAc]
                )
                self.simulations[observables.index('cfos_mRNA'), :, i] = (
                    Y[:, V.cfosmRNAc]
                )
                self.simulations[observables.index('cFos_Protein'), :, i] = (
                    (Y[:, V.pcFOSn] + Y[:, V.cFOSn]) * (x[C.Vn]/x[C.Vc])
                    + Y[:, V.cFOSc] + Y[:, V.pcFOSc]
                )
                self.simulations[observables.index('Phosphorylated_cFos'), :, i] = (
                    Y[:, V.pcFOSn] * (x[C.Vn]/x[C.Vc]) + Y[:, V.pcFOSc]
                )
    
    @staticmethod
    def _solveode(diffeq, y0, tspan, args, dt=1.0):
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
    standard_error = [None] * len(observables)

    t2 = [0, 300, 600, 900, 1800, 2700, 3600, 5400]  # (Unit: sec.)

    experiments[observables.index('Phosphorylated_MEKc')] = {
        'EGF': [0.000, 0.773, 0.439, 0.252, 0.130, 0.087, 0.080, 0.066], 
        'HRG': [0.000, 0.865, 1.000, 0.837, 0.884, 0.920, 0.875, 0.789], 
    }
    standard_error[observables.index('Phosphorylated_MEKc')] = {
        'EGF': [sd/np.sqrt(3) for sd in [0.000, 0.030, 0.048, 0.009, 0.009, 0.017, 0.012, 0.008]], 
        'HRG': [sd/np.sqrt(3) for sd in [0.000, 0.041, 0.000, 0.051, 0.058, 0.097, 0.157, 0.136]], 
    }

    experiments[observables.index('Phosphorylated_ERKc')] = {
        'EGF': [0.000, 0.867, 0.799, 0.494, 0.313, 0.266, 0.200, 0.194], 
        'HRG': [0.000, 0.848, 1.000, 0.971, 0.950, 0.812, 0.747, 0.595], 
    }
    standard_error[observables.index('Phosphorylated_ERKc')] = {
        'EGF': [sd/np.sqrt(3) for sd in [0.000, 0.137, 0.188, 0.126, 0.096, 0.087, 0.056, 0.012]], 
        'HRG': [sd/np.sqrt(3) for sd in [0.000, 0.120, 0.000, 0.037, 0.088, 0.019, 0.093, 0.075]], 
    }

    experiments[observables.index('Phosphorylated_RSKw')] = {
        'EGF': [0, 0.814, 0.812, 0.450, 0.151, 0.059, 0.038, 0.030], 
        'HRG': [0, 0.953, 1.000, 0.844, 0.935, 0.868, 0.779, 0.558], 
    }
    standard_error[observables.index('Phosphorylated_RSKw')] = {
        'EGF': [sd/np.sqrt(3) for sd in [0, 0.064, 0.194, 0.030, 0.027, 0.031, 0.043, 0.051]], 
        'HRG': [sd/np.sqrt(3) for sd in [0, 0.230, 0.118, 0.058, 0.041, 0.076, 0.090, 0.077]], 
    }

    experiments[observables.index('Phosphorylated_cFos')] = {
        'EGF': [0, 0.060, 0.109, 0.083, 0.068, 0.049, 0.027, 0.017], 
        'HRG': [0, 0.145, 0.177, 0.158, 0.598, 1.000, 0.852, 0.431], 
    }
    standard_error[observables.index('Phosphorylated_cFos')] = {
        'EGF': [sd/np.sqrt(3) for sd in [0, 0.003, 0.021, 0.013, 0.016, 0.007, 0.003, 0.002]], 
        'HRG': [sd/np.sqrt(3) for sd in [0, 0.010, 0.013, 0.001, 0.014, 0.000, 0.077, 0.047]], 
    }

    # --------------------------------------------------------------------------
    t3 = [0, 600, 1800, 3600, 5400]  # (Unit: sec.)

    experiments[observables.index('Phosphorylated_CREBw')] = {
        'EGF': [0, 0.446, 0.030, 0.000, 0.000], 
        'HRG': [0, 1.000, 0.668, 0.460, 0.340], 
    }
    standard_error[observables.index('Phosphorylated_CREBw')] = {
        'EGF': [sd/np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]], 
        'HRG': [sd/np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]], 
    }
    # --------------------------------------------------------------------------
    t4 = [0, 600, 1200, 1800, 2700, 3600, 5400]  # (Unit: sec.)

    experiments[observables.index('cfos_mRNA')] = {
        'EGF': [0, 0.181, 0.476, 0.518, 0.174, 0.026, 0.000], 
        'HRG': [0, 0.353, 0.861, 1.000, 0.637, 0.300, 0.059], 
    }
    standard_error[observables.index('cfos_mRNA')] = {
        'EGF': [sd/np.sqrt(3) for sd in [0.017, 0.004, 0.044, 0.004, 0.023, 0.007, 0.008]], 
        'HRG': [sd/np.sqrt(3) for sd in [0.017, 0.006, 0.065, 0.044, 0.087, 0.023, 0.001]], 
    }
    # --------------------------------------------------------------------------
    t5 = [0, 900, 1800, 2700, 3600, 5400]  # (Unit: sec.)

    experiments[observables.index('cFos_Protein')] = {
        'EGF': [0, 0.078, 0.216, 0.240, 0.320, 0.235], 
        'HRG': [0, 0.089, 0.552, 0.861, 1.000, 0.698], 
    }
    standard_error[observables.index('cFos_Protein')] = {
        'EGF': [sd/np.sqrt(3) for sd in [0, 0.036, 0.028, 0.056, 0.071, 0.048]], 
        'HRG': [sd/np.sqrt(3) for sd in [0, 0.021, 0.042, 0.063, 0.000, 0.047]], 
    }
    
    experiments[observables.index('dusp_mRNA')] = {
        'EGF': [0.000, 0.177, 0.331, 0.214, 0.177, 0.231], 
        'HRG': [0.000, 0.221, 0.750, 1.000, 0.960, 0.934], 
    }
    standard_error[observables.index('dusp_mRNA')] = {
        'EGF': [sd/np.sqrt(3) for sd in [0.033, 0.060, 0.061, 0.032, 0.068, 0.050]], 
        'HRG': [sd/np.sqrt(3) for sd in [0.027, 0.059, 0.094, 0.124, 0.113, 0.108]], 
    }

    def get_timepoint(self, obs_idx):
        if obs_idx in [
            observables.index('Phosphorylated_MEKc'),
            observables.index('Phosphorylated_ERKc'),
            observables.index('Phosphorylated_RSKw'),
            observables.index('Phosphorylated_cFos'),
        ]:
            exp_t = self.t2

        elif obs_idx == observables.index('Phosphorylated_CREBw'):
            exp_t = self.t3

        elif obs_idx == observables.index('cfos_mRNA'):
            exp_t = self.t4

        elif obs_idx in [
            observables.index('cFos_Protein'),
            observables.index('dusp_mRNA'),
        ]:
            exp_t = self.t5

        return list(map(int, exp_t))


class Visualization(object):
    def __init__(self):
        self.cm = plt.cm.get_cmap('tab20')

        self.timecourse_options = [
            {
                'divided_by' : 1,  # to convert time unit. (e.g. sec -> min)
                'xlim' : (),
                'xticks' : None,
                'xlabel': 'Time',
                'ylim' : (),
                'yticks' : None,
                'ylabel': observables[i].replace('__', '\n').replace('_', ' '),
                'cmap' : [self.cm.colors[j] for j in range(20)],
                'shape' : Line2D.filled_markers,
                'dont_show' : [],  # conditions you don't want to plot
            } for i, _ in enumerate(observables)]
        
        self.obs_samefig = {
            'observables' : [],
            'condition' : None,
            'xlim' : (),
            'xticks' : None,
            'xlabel': 'Time',
            'ylim' : (),
            'yticks' : None,
            'ylabel': '',
            'cmap' : [self.cm.colors[j] for j in range(20)],
            'shape' : Line2D.filled_markers,
        }
    
        self.sensitivity_options = {
            'figsize' : (12, 5),
            'width' : 0.3,
            'cmap' : [self.cm.colors[j] for j in range(20)],
        }

    def get_timecourse_options(self):
        '''
        for i, _ in enumerate(observables):
            self.timecourse_options[i]['divided_by'] = 60  # sec. -> min.
            self.timecourse_options[i]['xlim'] = (-5, 95)
            self.timecourse_options[i]['xticks'] = [0, 30, 60, 90]
            self.timecourse_options[i]['xlabel'] = 'Time (min)'
            self.timecourse_options[i]['ylim'] = (-0.1, 1.3)
            self.timecourse_options[i]['yticks'] = [0.0, 0.3, 0.6, 0.9, 1.2]
            self.timecourse_options[i]['cmap'] = ['mediumblue', 'red']
            self.timecourse_options[i]['shape'] = ['D', 's']
            self.timecourse_options[i]['dont_show'] = []

        self.timecourse_options[
            observables.index('Phosphorylated_MEKc')
        ]['ylabel'] = 'Phosphorylated MEK\n(cytoplasm)'

        self.timecourse_options[
            observables.index('Phosphorylated_ERKc')
        ]['ylabel'] = 'Phosphorylated ERK\n(cytoplasm)'

        self.timecourse_options[
            observables.index('Phosphorylated_RSKw')
        ]['ylabel'] = 'Phosphorylated RSK\n(whole cell)'

        self.timecourse_options[
            observables.index('Phosphorylated_CREBw')
        ]['ylabel'] = 'Phosphorylated CREB\n(whole cell)'

        self.timecourse_options[
            observables.index('dusp_mRNA')
        ]['ylabel'] = r'$\it{dusp}$'+' mRNA\nexpression'

        self.timecourse_options[
            observables.index('cfos_mRNA')
        ]['ylabel'] = r'$\it{c}$'+'-'+r'$\it{fos}$'+' mRNA\nexpression'

        self.timecourse_options[
            observables.index('cFos_Protein')
        ]['ylabel'] = 'c-Fos Protein\nexpression'

        self.timecourse_options[
            observables.index('Phosphorylated_cFos')
        ]['ylabel'] = 'Phosphorylated c-Fos\nProtein expression'
        '''
        return self.timecourse_options

    def multiplot_observables(self):
        ''' Example
        self.obs_samefig['observables'] = [
            'Phosphorylated_ERKc',
            'Phosphorylated_CREBw',
            'Phosphorylated_cFos',
            'dusp_mRNA'
        ]
        self.obs_samefig['condition'] = 'EGF'
        self.obs_samefig['xlim'] = (-5, 95)
        self.obs_samefig['xticks'] = [0, 30, 60, 90]
        self.obs_samefig['xlabel'] = 'Time (min)'
        self.obs_samefig['ylim'] = (-0.1, 1.3)
        self.obs_samefig['yticks'] = [0.0, 0.3, 0.6, 0.9, 1.2]
        self.obs_samefig['ylabel'] = 'Intensity (a.u.)'
        '''
        return self.obs_samefig
    
    @staticmethod
    def set_timecourse_rcParams():
        """ figure/simulation
        """
        plt.rcParams['font.size'] = 20
        plt.rcParams['axes.linewidth'] = 1.5
        plt.rcParams['xtick.major.width'] = 1.5
        plt.rcParams['ytick.major.width'] = 1.5
        plt.rcParams['lines.linewidth'] = 1.8
        plt.rcParams['lines.markersize'] = 12
        # plt.rcParams['font.family'] = 'Arial'
        # plt.rcParams['mathtext.fontset'] = 'custom'
        # plt.rcParams['mathtext.it'] = 'Arial:italic'

    @staticmethod
    def set_param_range_rcParams():
        """ figure/param_range
        """
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def set_sensitivity_rcParams():
        """ figure/sensitivity
        """
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def convert_species_name(name):
        """ figure/sensitivity/initial_condition
        - Sensitivity for species with nonzero initial conditions
        """
        '''
        if name == 'ERKc':
            return 'ERK (cytoplasm)'
        elif name == 'RSKc':
            return 'RSK (cytoplasm)'
        elif name == 'CREBn':
            return 'CREB (nucleus)'
        elif name == 'Elk1n':
            return 'Elk1 (nucleus)'
        '''
        return name