import os
import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from biomass.model.name2idx import variables as V
from biomass.model.initial_condition import initial_values
from .sensitivity import analyze_sensitivity

os.makedirs('./figure/sensitivity/nonzero_init', exist_ok=True)

width = 0.3

nonzero_idx = []
y0 = initial_values()
for i,val in enumerate(y0):
    if val != 0.0:
        nonzero_idx.append(i)
if len(nonzero_idx) == 0:
    print('No nonzero initial values')
    sys.exit()

if not os.path.isfile('sensitivities_npy/nonzero_init/s_cFosmRNA.npy') or \
    not os.path.isfile('sensitivities_npy/nonzero_init/s_PcFos.npy'):
        
    os.makedirs('./sensitivities_npy/nonzero_init', exist_ok=True)
    (s_cFosmRNA, s_PcFos) = analyze_sensitivity(nonzero_idx)
    np.save('sensitivities_npy/nonzero_init/s_cFosmRNA',s_cFosmRNA)
    np.save('sensitivities_npy/nonzero_init/s_PcFos',s_PcFos)
else:
    s_cFosmRNA = np.load('sensitivities_npy/nonzero_init/s_cFosmRNA.npy')
    s_PcFos = np.load('sensitivities_npy/nonzero_init/s_PcFos.npy')
            

def sensitivity_barplot():
    plt.figure(figsize=(9,5))
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['axes.linewidth'] = 1

    average = np.nanmean(s_cFosmRNA,axis=0)
    stdev = np.nanstd(s_cFosmRNA,axis=0,ddof=1)

    plt.hlines([0],-width,len(nonzero_idx)-1-width,'k',lw=1)
    plt.bar(
        np.arange(len(nonzero_idx)),average[:-1,0],
        yerr = stdev[:-1,0],ecolor='mediumblue',capsize=2,
        width=width,color='mediumblue',align='center',label='EGF'
    )
    plt.bar(
        np.arange(len(nonzero_idx))+width,average[:-1,1],
        yerr = stdev[:-1,1],ecolor='red',capsize=2,
        width=width,color='red',align='center',label='HRG'
    )

    plt.xticks(np.arange(len(nonzero_idx))+width/2,[V.var_names[i] for i in nonzero_idx],rotation=30)
    plt.ylabel('Control coefficients on\nduration ('+r'$\it{c}$'+'-'+r'$\it{fos}$'+' mRNA)')
    plt.xlim(-width,len(nonzero_idx)-1-width)
    plt.ylim(-0.2,1.2)
    plt.legend(loc='upper left',frameon=False)
    plt.savefig('figure/sensitivity/nonzero_init/cFosmRNA.pdf',bbox_inches='tight')
    plt.close()

    # ==========================================================================
    
    plt.figure(figsize=(9,5))
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['axes.linewidth'] = 1
    
    average = np.nanmean(s_PcFos,axis=0)
    stdev = np.nanstd(s_PcFos,axis=0,ddof=1)

    plt.hlines([0],-width,len(nonzero_idx)-1-width,'k',lw=1)
    plt.bar(
        np.arange(len(nonzero_idx)),average[:-1,0],
        yerr = stdev[:-1,0],ecolor='mediumblue',capsize=2,
        width=width,color='mediumblue',align='center',label='EGF'
    )
    plt.bar(
        np.arange(len(nonzero_idx))+width,average[:-1,1],
        yerr = stdev[:-1,1],ecolor='red',capsize=2,
        width=width,color='red',align='center',label='HRG'
    )

    plt.xticks(np.arange(len(nonzero_idx))+width/2,[V.var_names[i] for i in nonzero_idx],rotation=30)
    plt.ylabel('Control coefficients on\nintegrated response (pc-Fos)')
    plt.xlim(-width,len(nonzero_idx)-1-width)
    plt.legend(loc='upper left',frameon=False)

    plt.savefig('figure/sensitivity/nonzero_init/PcFos.pdf',bbox_inches='tight')
    plt.close()


def sensitivity_heatmap():
    if len(nonzero_idx) < 2:
        pass
    else:
        # e.g. Sensitivity coefficients on duration (c-fos mRNA, HRG-induced)
        sensitivity_matrix = s_cFosmRNA[:,:,1]
        
        # Normalize from -1 to 1
        nanidx =[]
        for i in range(sensitivity_matrix.shape[0]):
            if any(np.isnan(sensitivity_matrix[i,:])):
                nanidx.append(i)
            sensitivity_matrix[i,:] = \
                sensitivity_matrix[i,:]/np.nanmax(np.abs(sensitivity_matrix[i,:]))
        sensitivity_matrix = np.delete(sensitivity_matrix,nanidx,axis=0)
        
        plt.rcParams['font.size'] = 8
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = 'Arial:italic'
        plt.rcParams['axes.linewidth'] = 1
        
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
        
        
        plt.suptitle('Control coefficients on duration ('+r'$\it{c}$'+'-'+r'$\it{fos}$'+' mRNA)', fontsize=24)
        plt.savefig('figure/sensitivity/nonzero_init/h_cFosmRNA_hrg.pdf',bbox_inches='tight')
        plt.close()
        
        # ==========================================================================
        
        # e.g. Sensitivity coefficients on integrated response (pc-Fos, EGF-induced)
        sensitivity_matrix = s_PcFos[:,:,0]
        
        # Normalize from -1 to 1
        nanidx =[]
        for i in range(sensitivity_matrix.shape[0]):
            if any(np.isnan(sensitivity_matrix[i,:])):
                nanidx.append(i)
            sensitivity_matrix[i,:] = \
                sensitivity_matrix[i,:]/np.nanmax(np.abs(sensitivity_matrix[i,:]))
        sensitivity_matrix = np.delete(sensitivity_matrix,nanidx,axis=0)
        
        plt.rcParams['font.size'] = 8
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = 'Arial:italic'
        plt.rcParams['axes.linewidth'] = 1

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
        
        plt.savefig('figure/sensitivity/nonzero_init/h_PcFos_egf.pdf',bbox_inches='tight')
        plt.close()