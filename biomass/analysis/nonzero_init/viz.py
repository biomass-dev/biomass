import os
import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from biomass.model.name2idx import variables as V
from biomass.model.initial_condition import initial_values
from biomass.observable import observables, NumericalSimulation
from .sensitivity import analyze_sensitivity

os.makedirs('./figure/sensitivity/nonzero_init/heatmap', exist_ok=True)

sim = NumericalSimulation()

width = 0.3

nonzero_idx = []
y0 = initial_values()
for i,val in enumerate(y0):
    if val != 0.0:
        nonzero_idx.append(i)
if len(nonzero_idx) == 0:
    print('No nonzero initial values')
    sys.exit()

def run_analysis(metric):
    if not os.path.isfile(
        'sensitivities_npy/nonzero_init/%s/sensitivity_coefficients.npy'%(metric)):
        os.makedirs('./sensitivities_npy/nonzero_init/%s'%(metric),exist_ok=True)
        sensitivity_coefficients = analyze_sensitivity(metric,nonzero_idx)
        np.save(
            'sensitivities_npy/nonzero_init/%s/sensitivity_coefficients'%(metric),
            sensitivity_coefficients
        )
    else:
        sensitivity_coefficients = np.load(
            'sensitivities_npy/nonzero_init/%s/sensitivity_coefficients.npy'%(metric)
        )
    return sensitivity_coefficients


def sensitivity_barplot(metric):
    sensitivity_coefficients = run_analysis(metric)
    
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rcParams['xtick.major.width'] = 1.2
    plt.rcParams['ytick.major.width'] = 1.2

    color = ['mediumblue','red']
    for k,obs_name in enumerate(observables):
        plt.figure(figsize=(9,5))
        plt.hlines([0],-width,len(nonzero_idx)-1-width,'k',lw=1)
        for l,condition in enumerate(sim.conditions):
            sensitivity_matrix = sensitivity_coefficients[:,:,k,l]
            nan_idx = []
            for m in range(sensitivity_matrix.shape[0]):
                if any(np.isnan(sensitivity_matrix[m,:])):
                    nan_idx.append(m)
            sensitivity_matrix = np.delete(
                sensitivity_matrix, nan_idx, axis=0
            )
            if sensitivity_matrix.size != 0:
                average = np.mean(sensitivity_matrix,axis=0)
                stdev = np.std(sensitivity_matrix,axis=0,ddof=1)
                plt.bar(
                    np.arange(len(nonzero_idx))+l*width, average, yerr = stdev,
                    ecolor=color[l],capsize=2,width=width,color=color[l],
                    align='center',label=condition
                )
        plt.xticks(
            np.arange(len(nonzero_idx))+width/2,[V.var_names[i] for i in nonzero_idx],
            rotation=30
        )
        plt.ylabel(
            'Control coefficients on\n'+metric+' ('+obs_name.replace('_',' ')+')'
        )
        plt.xlim(-width,len(nonzero_idx)-1-width)
        plt.legend(loc='upper left',frameon=False)
        plt.savefig(
            'figure/sensitivity/nonzero_init/{0}_{1}.pdf'.format(obs_name,metric),
            bbox_inches='tight'
        )
        plt.close()


def sensitivity_heatmap(metric):
    if len(nonzero_idx) < 2:
        pass
    else:
        sensitivity_coefficients = run_analysis(metric)
        
        plt.rcParams['font.size'] = 8
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = 'Arial:italic'
                
        for k,obs_name in enumerate(observables):
            for l,condition in enumerate(sim.conditions):
                sensitivity_matrix = sensitivity_coefficients[:,:,k,l]
                # Normalize from -1 to 1
                nan_idx = []
                for i in range(sensitivity_matrix.shape[0]):
                    if any(np.isnan(sensitivity_matrix[i,:])):
                        nan_idx.append(i)
                    if np.nanmax(np.abs(sensitivity_matrix[i,:])) == 0.0:
                        sensitivity_matrix[i,:] = \
                            np.zeros(sensitivity_matrix.shape[1])
                    else:
                        sensitivity_matrix[i,:] = \
                            sensitivity_matrix[i,:]/np.nanmax(np.abs(sensitivity_matrix[i,:]))
                sensitivity_matrix = np.delete(
                    sensitivity_matrix, nan_idx, axis=0
                )
                if sensitivity_matrix.size != 0:
                    sns.clustermap(
                        sensitivity_matrix,
                        center=0,
                        method='ward',
                        cmap='RdBu_r',
                        linewidth=.5,
                        col_cluster=False,
                        figsize = (16,8),
                        xticklabels=[V.var_names[i] for i in nonzero_idx],
                        yticklabels=[],
                        cbar_kws={"ticks":[-1,0,1]}
                    )
                    plt.savefig(
                        'figure/sensitivity/nonzero_init/heatmap/{0}_{1}_{2}.pdf'.format(condition,obs_name,metric),
                        bbox_inches='tight'
                    )
                    plt.close()