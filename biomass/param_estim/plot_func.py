import numpy as np
from matplotlib import pyplot as plt

from biomass.observable import observable_names, num_observables, ExperimentalData

def timecourse(sim,n_file,viz_type,show_all,stdev,simulations_all):

    exp = ExperimentalData()

    for j,title in enumerate(observable_names):

        plt.figure(figsize=(4,3))
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.rcParams['font.size'] = 12
        """
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = 'Arial:italic'
        """
        plt.rcParams['axes.linewidth'] = 1
        plt.rcParams['lines.linewidth'] = 2
        plt.rcParams['lines.markersize'] = 10

        if show_all:
            for i in range(n_file):
                plt.plot(sim.t,simulations_all[j,i,:,0]/np.max(simulations_all[j,i,:,:]),'b',alpha=0.05)
                plt.plot(sim.t,simulations_all[j,i,:,1]/np.max(simulations_all[j,i,:,:]),'r',alpha=0.05)

        if not viz_type == 'average':
            plt.plot(sim.t,sim.simulations[j,:,0]/np.max(sim.simulations[j]),'b')
            plt.plot(sim.t,sim.simulations[j,:,1]/np.max(sim.simulations[j]),'r')
        else:
            normalized = np.empty((num_observables,n_file,len(sim.tspan),sim.condition))
            for i in range(n_file):
                normalized[j,i,:,0] = simulations_all[j,i,:,0]/np.max(simulations_all[j,i,:,:])
                normalized[j,i,:,1] = simulations_all[j,i,:,1]/np.max(simulations_all[j,i,:,:])
            plt.plot(sim.t,np.nanmean(normalized[j,:,:,0],axis=0),'b')
            plt.plot(sim.t,np.nanmean(normalized[j,:,:,1],axis=0),'r')
            if stdev:
                mean_egf = np.nanmean(normalized[j,:,:,0],axis=0)
                yerr_egf = [np.nanstd(normalized[j,:,i,0],ddof=1) for i,_ in enumerate(sim.t)]
                plt.fill_between(
                    sim.t, mean_egf - yerr_egf, mean_egf + yerr_egf,
                    lw=0,color='b',alpha=0.1
                )
                mean_hrg = np.nanmean(normalized[j,:,:,1],axis=0)
                yerr_hrg = [np.nanstd(normalized[j,:,i,1],ddof=1) for i,_ in enumerate(sim.t)]
                plt.fill_between(
                    sim.t, mean_hrg - yerr_hrg, mean_hrg + yerr_hrg,
                    lw=0,color='r',alpha=0.1
                )

        if exp.experiments[j] is not None:
            exp_t = exp.get_timepoint(j)
            plt.plot(exp_t/60.,exp.experiments[j]['EGF'],'D',markerfacecolor='None',markeredgecolor='b',clip_on=False)
            plt.plot(exp_t/60.,exp.experiments[j]['HRG'],'s',markerfacecolor='None',markeredgecolor='r',clip_on=False)

        plt.xlim(0,90)
        plt.xticks([0,30,60,90])
        plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
        plt.ylim(0,1.2)
        plt.xlabel('Time (min)')
        plt.title(title)

        plt.savefig('./figure/simulation_{0}_{1}.pdf'.
                    format(viz_type,title),bbox_inches='tight')
        plt.close()