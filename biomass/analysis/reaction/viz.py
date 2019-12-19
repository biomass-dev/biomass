import os
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from biomass.observable import observables,NumericalSimulation
from .sensitivity import analyze_sensitivity
from .reaction import *

os.makedirs('./figure/sensitivity/reaction/heatmap', exist_ok=True)

sim = NumericalSimulation()
width = 0.3

def run_analysis(metric):
    if not os.path.isfile(
        'sensitivities_npy/reaction/%s/sensitivity_coefficients.npy'%(metric)):
        os.makedirs('./sensitivities_npy/reaction/%s'%(metric),exist_ok=True)
        sensitivity_coefficients = analyze_sensitivity(metric,num_reaction)
        np.save(
            'sensitivities_npy/reaction/%s/sensitivity_coefficients'%(metric),
            sensitivity_coefficients
        )
    else:
        sensitivity_coefficients = np.load(
            'sensitivities_npy/reaction/%s/sensitivity_coefficients.npy'%(metric)
        )
    return sensitivity_coefficients


def get_sort_idx():
    reaction_module = get_reaction_module()
    sort_idx = [0]*num_reaction
    left_end = 0
    for i,ith_module in enumerate(reaction_module):
        for j,k in enumerate(ith_module):
            if i != 0 and j == 0:
                left_end += len(reaction_module[i-1])
            sort_idx[left_end+j] = k
    return sort_idx


def get_reaction_number(sort_idx):
    
    return [str(i) for i in sort_idx]


def draw_vertical_span(width):
    reaction_module = get_reaction_module()
    left_end = 0
    for i,ith_module in enumerate(reaction_module):
        if i%2 == 0:
            plt.axvspan(
                left_end-width,
                left_end+len(ith_module)-width,
                facecolor='k',alpha=0.1
            )
        left_end += len(ith_module)


def sensitivity_barplot(metric):
    sensitivity_coefficients = run_analysis(metric)
    reaction_module = get_reaction_module()
    sort_idx = get_sort_idx()
    reaction_number = get_reaction_number(sort_idx)
    
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rcParams['xtick.major.width'] = 1.2
    plt.rcParams['ytick.major.width'] = 1.2
    
    color = ['mediumblue','red']
    for obs_name in observables:
        k = observables.index(obs_name)
        average = np.nanmean(sensitivity_coefficients[:,:,k,:],axis=0)
        stdev = np.nanstd(sensitivity_coefficients[:,:,k,:],axis=0,ddof=1)

        plt.figure(figsize=(12,5))
        draw_vertical_span(width)
        plt.hlines([0],-width,num_reaction-1-width,'k',lw=1)
        for l,condition in enumerate(sim.conditions):
            plt.bar(
                np.arange(num_reaction)+l*width,average[sort_idx,l],
                yerr = stdev[sort_idx,l],
                ecolor=color[l],capsize=2,width=width,color=color[l],
                align='center',label=condition
            )
        distance = np.max(average)*0.05
        for i,j in enumerate(sort_idx):
            if j != 0:
                xp = i + width/2
                yp = average[j,np.argmax(np.abs(average[j,:]))]
                yerr = stdev[j,np.argmax(stdev[j,:])]
                if yp > 0:
                    plt.text(
                        xp,yp+yerr+distance,reaction_number[i],
                        ha='center', va='bottom', fontsize=10, rotation=90
                    )
                else:
                    plt.text(
                        xp,yp-yerr-distance,reaction_number[i],
                        ha='center', va='top', fontsize=10, rotation=90
                    )
        plt.xticks([])
        plt.ylabel(
            'Control coefficients on\n'+metric+' ('+obs_name.replace('_',' ')+')'
        )
        plt.xlim(-width,num_reaction-1-width)
        #plt.ylim(-1.2,0.6)
        #plt.yticks([-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6])
        plt.legend(loc='lower right',frameon=False)
        plt.savefig(
            'figure/sensitivity/reaction/{0}_{1}.pdf'.format(obs_name,metric),
            bbox_inches='tight'
        )
        plt.close()

    
def sensitivity_heatmap(metric):
    sensitivity_coefficients = run_analysis(metric)
    reaction_module = get_reaction_module()
    sort_idx = get_sort_idx()
    reaction_number = get_reaction_number(sort_idx)
    
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rcParams['xtick.major.width'] = 1.2
    plt.rcParams['ytick.major.width'] = 1.2
    
    for obs_name in observables:
        k = observables.index(obs_name)
        for l,condition in enumerate(sim.conditions):
            sensitivity_matrix = sensitivity_coefficients[:,sort_idx[:-1],k,l]
            # Normalize from -1 to 1
            nanidx =[]
            for i in range(sensitivity_matrix.shape[0]):
                if any(np.isnan(sensitivity_matrix[i,:])):
                    nanidx.append(i)
                sensitivity_matrix[i,:] = \
                    sensitivity_matrix[i,:]/np.nanmax(np.abs(sensitivity_matrix[i,:]))
            sensitivity_matrix = np.delete(sensitivity_matrix,nanidx,axis=0)
            
            sns.clustermap(
                sensitivity_matrix,
                center=0,
                method='ward',
                cmap='RdBu_r',
                linewidth=.5,
                col_cluster=False,
                figsize = (16,8),
                xticklabels=[reaction_number[i] for i in range(num_reaction-1)],
                yticklabels=[],
                cbar_kws={"ticks":[-1,0,1]}
            )
            plt.savefig(
                'figure/sensitivity/reaction/heatmap/{0}_{1}_{2}.pdf'.format(condition,obs_name,metric),
                bbox_inches='tight'
            )
            plt.close()