import os
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from biomass import model
from biomass.observable import species, NumericalSimulation
from biomass.param_estim import plot_func
from biomass.param_estim.search_parameter import search_parameter_index

def simulate_all(viz_type,show_all,stdev):
    """Simulate ODE model with estimated parameter values.
    
    Parameters
    ----------
    viz_type : str
        - 'average': The average of simulation results with parameter sets in "out/"
        - 'best': The best simulation result in "out/", simulation with "best_fit_param"
        - 'original': Simulation with the default parameters and initial values defined in "biomass/model/"
        - 'n(=1,2,...)': Use the parameter set in "out/n/"
    show_all : bool
        Whether to show all simulation results
    stdev: bool
        If True, the standard deviation of simulated values will be shown
        (only when viz_type == 'average')
    """
    if not viz_type in ['best','average','original']:
        try:
            int(viz_type)
        except ValueError:
            print("Error: viz_type ∈ {'best','average','original',int(1~n_fit_param)}")

    x = model.f_params()
    y0 = model.initial_values()
    sim = NumericalSimulation()

    n_file = 0
    if viz_type == 'original':
        pass
    else:
        if os.path.isdir('./out'):
            fit_param_files = os.listdir('./out')
            for file in fit_param_files:
                if re.match(r'\d',file):
                    n_file += 1

    simulations_all = np.full((len(species),n_file,len(sim.t),len(sim.conditions)),np.nan)
    
    if n_file > 0:
        if n_file == 1 and viz_type == 'average':
            viz_type = 'best'
        for i in range(n_file):
            (sim,successful) = validate(i+1,x,y0)
            if successful:
                for j,_ in enumerate(species):
                    simulations_all[j,i,:,:] = sim.simulations[j,:,:]

        best_fitness_all = np.empty(n_file)
        for i in range(n_file):
            if os.path.isfile('./out/%d/best_fitness.npy'%(i+1)):
                best_fitness_all[i] = np.load('./out/%d/best_fitness.npy'%(i+1))
            else:
                best_fitness_all[i] = np.inf

        best_paramset = np.argmin(best_fitness_all) + 1
        write_best_fit_param(best_paramset,x,y0)

        if viz_type == 'average':
            pass
        elif viz_type == 'best':
            sim = validate(int(best_paramset),x,y0)[0]
        elif int(viz_type) <= n_file:
            sim = validate(int(viz_type),x,y0)[0]
        else:
            raise ValueError(
                '%d is larger than n_fit_param(%d)'%(int(viz_type),n_file)
            )

        if n_file >= 2:
            save_param_range(n_file,x,y0,portrait=True)

    else:
        if sim.simulate(x,y0) is not None:
            print('Simulation failed.')

    plot_func.timecourse(sim,n_file,viz_type,show_all,stdev,simulations_all)


def update_param(paramset,x,y0):
    search_idx = search_parameter_index()

    if os.path.isfile('./out/%d/generation.npy'%(paramset)):
        generation = np.load('./out/%d/generation.npy'%(paramset))
        best_indiv = np.load('./out/%d/fit_param%d.npy'%(paramset,int(generation)))

        for i,j in enumerate(search_idx[0]):
            x[j] = best_indiv[i]
        for i,j in enumerate(search_idx[1]):
            y0[j] = best_indiv[i+len(search_idx[0])]
    else:
        pass
    return x, y0


def validate(nth_paramset,x,y0):
    # -------------------------------------------------------------------------
    # Validates the dynamical viability of a set of estimated parameter values.
    # -------------------------------------------------------------------------
    sim = NumericalSimulation()
    
    (x,y0) = update_param(nth_paramset,x,y0)

    if sim.simulate(x,y0) is None:
        return sim,True
    else:
        print('Simulation failed.\nparameter_set #%d'%(nth_paramset))
        return sim,False


def write_best_fit_param(best_paramset,x,y0):

    (x,y0) = update_param(best_paramset,x,y0)
    
    with open('./out/best_fit_param.txt', mode='w') as f:
        f.write('# param set: %d\n'%(best_paramset))
        f.write('\n### Param. const\n')
        for i in range(model.C.len_f_params):
            f.write('x[C.%s] = %8.3e\n'%(model.C.param_names[i],x[i]))
        f.write('\n### Non-zero initial conditions\n')
        for i in range(model.V.len_f_vars):
            if y0[i] != 0:
                f.write('y0[V.%s] = %8.3e\n'%(model.V.var_names[i],y0[i]))


def save_param_range(n_file,x,y0,portrait):
    search_idx = search_parameter_index()
    search_param_matrix = np.empty((n_file,len(search_idx[0])+len(search_idx[1])))

    for nth_paramset in range(1,n_file+1):
        if os.path.isfile('./out/%d/generation.npy'%(nth_paramset)):
            generation = np.load('./out/%d/generation.npy'%(nth_paramset))
            best_indiv = np.load('./out/%d/fit_param%d.npy'%(nth_paramset,int(generation)))
        else:
            best_indiv = np.empty(len(search_idx[0])+len(search_idx[1]))
            for i,j in enumerate(search_idx[0]):
                best_indiv[i] = x[j]
            for i,j in enumerate(search_idx[1]):
                best_indiv[i+len(search_idx[0])] = y0[j]

        search_param_matrix[nth_paramset-1,:] = best_indiv

    # --------------------------------------------------------------------------
    # seaborn.boxenplot
    if portrait:
        if len(search_idx[0]) > 0:
            fig = plt.figure(figsize=(8,len(search_idx[0])/2.5))
            plt.gca().spines['right'].set_visible(False)
            plt.gca().spines['top'].set_visible(False)
            plt.gca().yaxis.set_ticks_position('left')
            plt.gca().xaxis.set_ticks_position('bottom')

            ax = sns.boxenplot(
                data=search_param_matrix[:,:len(search_idx[0])],
                orient='h',
                linewidth=0.5,
                palette='Set2'
            )
            ax.set_xlabel('Parameter value')
            ax.set_ylabel('')
            ax.set_yticklabels([model.C.param_names[i] for i in search_idx[0]])
            ax.set_xscale('log')

            plt.savefig('./figure/param_range.pdf',bbox_inches='tight')
            plt.close(fig)
        if len(search_idx[1]) > 0:
            fig = plt.figure(figsize=(8,len(search_idx[1])/2.5))
            plt.gca().spines['right'].set_visible(False)
            plt.gca().spines['top'].set_visible(False)
            plt.gca().yaxis.set_ticks_position('left')
            plt.gca().xaxis.set_ticks_position('bottom')

            ax = sns.boxenplot(
                data=search_param_matrix[:,len(search_idx[0]):],
                orient='h',
                linewidth=0.5,
                palette='Set2'
            )
            ax.set_xlabel('Initial value')
            ax.set_ylabel('')
            ax.set_yticklabels([model.V.var_names[i] for i in search_idx[1]])
            ax.set_xscale('log')

            plt.savefig('./figure/initial_value_range.pdf',bbox_inches='tight')
            plt.close(fig)
    else:
        if len(search_idx[0]) > 0:
            fig = plt.figure(figsize=(len(search_idx[0])/2.2,6))
            plt.gca().spines['right'].set_visible(False)
            plt.gca().spines['top'].set_visible(False)
            plt.gca().yaxis.set_ticks_position('left')
            plt.gca().xaxis.set_ticks_position('bottom')

            ax = sns.boxenplot(
                data=search_param_matrix[:,:len(search_idx[0])],
                linewidth=0.5,
                palette='Set2'
            )
            ax.set_xlabel('')
            ax.set_xticklabels([model.C.param_names[i] for i in search_idx[0]],rotation=45)
            ax.set_ylabel('Parameter value')
            ax.set_yscale('log')

            plt.savefig('./figure/param_range_h.pdf',bbox_inches='tight')
            plt.close(fig)
        if len(search_idx[1]) > 0:
            fig = plt.figure(figsize=(len(search_idx[1])/2.2,6))
            plt.gca().spines['right'].set_visible(False)
            plt.gca().spines['top'].set_visible(False)
            plt.gca().yaxis.set_ticks_position('left')
            plt.gca().xaxis.set_ticks_position('bottom')

            ax = sns.boxenplot(
                data=search_param_matrix[:,len(search_idx[0]):],
                linewidth=0.5,
                palette='Set2'
            )
            ax.set_xlabel('')
            ax.set_xticklabels([model.V.var_names[i] for i in search_idx[1]],rotation=45)
            ax.set_ylabel('Initial value')
            ax.set_yscale('log')

            plt.savefig('./figure/initail_value_range_h.pdf',bbox_inches='tight')
            plt.close(fig)