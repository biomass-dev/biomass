import os
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from .sensitivity import analyze_sensitivity


def get_reaction_module():

    n_module = 16
    reaction_module = [None]*n_module

    # ERK_activation
    reaction_module[0] = [i for i in range(1,7)]

    # ERK_dephosphorylation_by_DUSP
    reaction_module[1] = [i for i in range(47,57)]

    # ERK_transport
    reaction_module[2] = [i for i in range(7,10)]

    # RSK_activation
    reaction_module[3] = [24,25]

    # RSK_transport
    reaction_module[4] = [26]

    # Elk1_activation
    reaction_module[5] = [29,30]

    # CREB_activation
    reaction_module[6] = [27,28]

    # dusp_production_etc
    reaction_module[7] = [i for i in range(10,14)]

    # DUSP_transport
    reaction_module[8] = [18,19]

    # DUSP_stabilization
    reaction_module[9] = [14,15,20,21]

    # DUSP_degradation
    reaction_module[10] = [16,17,22,23]

    # cfos_production_etc
    reaction_module[11] = [i for i in range(31,35)]

    # cFos_transport
    reaction_module[12] = [40,41]

    # cFos_stabilization
    reaction_module[13] = [35,36,37,42,43,44]

    # cFos_degradation
    reaction_module[14] = [38,39,45,46]
    
    # Feedback_from_F
    reaction_module[15] = [i for i in range(57,64)]

    return reaction_module


num_reaction = 64  # Num. of Rate Equations
width = 0.3

if not os.path.isfile('sensitivities_npy/s_cFosmRNA.npy') or not os.path.isfile('sensitivities_npy/s_PcFos.npy'):
    os.makedirs('./sensitivities_npy', exist_ok=True)
    (s_cFosmRNA, s_PcFos) = analyze_sensitivity(num_reaction)
    np.save('sensitivities_npy/s_cFosmRNA',s_cFosmRNA)
    np.save('sensitivities_npy/s_PcFos',s_PcFos)
else:
    s_cFosmRNA = np.load('sensitivities_npy/s_cFosmRNA.npy')
    s_PcFos = np.load('sensitivities_npy/s_PcFos.npy')
    
reaction_module = get_reaction_module()

sort_idx = [0]*num_reaction
left_end = 0
for i,ith_module in enumerate(reaction_module):
    for j,k in enumerate(ith_module):
        if i != 0 and j == 0:
            left_end += len(reaction_module[i-1])
        sort_idx[left_end+j] = k

reaction_number = [str(i) for i in sort_idx]


def draw_vertical_span(reaction_module,num_reaction,width):
    left_end = 0
    for i,ith_module in enumerate(reaction_module):
        if i%2 == 0:
            plt.axvspan(
                left_end-width,
                left_end+len(ith_module)-width,
                facecolor='k',alpha=0.1
            )
        left_end += len(ith_module)
        

def sensitivity_barplot():
    plt.figure(figsize=(12,5))
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['axes.linewidth'] = 1

    average = np.nanmean(s_cFosmRNA,axis=0)
    stdev = np.nanstd(s_cFosmRNA,axis=0,ddof=1)

    draw_vertical_span(reaction_module,num_reaction,width)
    plt.hlines([0],-width,num_reaction-1-width,'k',lw=1)
    plt.bar(
        np.arange(num_reaction),average[sort_idx,0],
        yerr = stdev[sort_idx,0],ecolor='mediumblue',capsize=2,
        width=width,color='mediumblue',align='center',label='EGF'
    )
    plt.bar(
        np.arange(num_reaction)+width,average[sort_idx,1],
        yerr = stdev[sort_idx,1],ecolor='red',capsize=2,
        width=width,color='red',align='center',label='HRG'
    )

    for i,j in enumerate(sort_idx):
        if j != 0:
            xp = i + width/2
            yp = average[j,np.argmax(np.abs(average[j,:]))]
            yerr = stdev[j,np.argmax(stdev[j,:])]
            if yp > 0:
                plt.text(
                    xp,yp+yerr+0.05,reaction_number[i],
                    ha='center', va='bottom', fontsize=10, rotation=90
                )
            else:
                plt.text(
                    xp,yp-yerr-0.05,reaction_number[i],
                    ha='center', va='top', fontsize=10, rotation=90
                )

    plt.xticks([])
    plt.ylabel('Control coefficients on\nduration ('+r'$\it{c}$'+'-'+r'$\it{fos}$'+' mRNA)')
    plt.xlim(-width,num_reaction-1-width)
    plt.ylim(-1.2,0.6)
    plt.yticks([-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6])
    plt.legend(loc='lower right',frameon=False)
    plt.savefig('figure/sensitivity_cFosmRNA.pdf',bbox_inches='tight')
    plt.close()

    # ==========================================================================
    
    plt.figure(figsize=(12,5))
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['axes.linewidth'] = 1
    
    average = np.nanmean(s_PcFos,axis=0)
    stdev = np.nanstd(s_PcFos,axis=0,ddof=1)

    draw_vertical_span(reaction_module,num_reaction,width)
    plt.hlines([0],-width,num_reaction-1-width,'k',lw=1)
    plt.bar(
        np.arange(num_reaction),average[sort_idx,0],
        yerr = stdev[sort_idx,0],ecolor='mediumblue',capsize=2,
        width=width,color='mediumblue',align='center',label='EGF'
    )
    plt.bar(
        np.arange(num_reaction)+width,average[sort_idx,1],
        yerr = stdev[sort_idx,1],ecolor='red',capsize=2,
        width=width,color='red',align='center',label='HRG'
    )

    for i,j in enumerate(sort_idx):
        if j != 0:
            xp = i + width/2
            yp = average[j,np.argmax(np.abs(average[j,:]))]
            yerr = stdev[j,np.argmax(stdev[j,:])]
            if yp > 0:
                plt.text(
                    xp,yp+yerr+0.2,reaction_number[i],
                    ha='center', va='bottom', fontsize=10, rotation=90
                )
            else:
                plt.text(
                    xp,yp-yerr-0.2,reaction_number[i],
                    ha='center', va='top', fontsize=10, rotation=90
                )

    plt.xticks([])
    plt.ylabel('Control coefficients on\nintegrated response (pc-Fos)')
    plt.xlim(-width,num_reaction-1-width)
    plt.ylim(-7,4.5)
    plt.legend(loc='lower right',frameon=False)

    plt.savefig('figure/sensitivity_PcFos.pdf',bbox_inches='tight')
    
    
def sensitivity_heatmap():
    # e.g. Sensitivity coefficients on duration (c-fos mRNA, HRG-induced)
    sensitivity_matrix = s_cFosmRNA[:,sort_idx[:-1],1]
    
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
        xticklabels=[reaction_number[i] for i in range(num_reaction-1)],
        yticklabels=[],
        cbar_kws={"ticks":[-1,0,1]}
    )
    
    plt.suptitle(
        'Normalized sensitivity coefficients on duration ('+r'$\it{c}$'+'-'+r'$\it{fos}$'+' mRNA)',
        fontsize=24
    )
    plt.savefig('figure/sensitivity_h_cFosmRNA_hrg.pdf',bbox_inches='tight')
    plt.close()
    
    # ==========================================================================
    
    # e.g. Sensitivity coefficients on integrated response (pc-Fos, EGF-induced)
    sensitivity_matrix = s_PcFos[:,sort_idx[:-1],0]
    
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
        xticklabels=[reaction_number[i] for i in range(num_reaction-1)],
        yticklabels=[],
        cbar_kws={"ticks":[-1,0,1]}
    )
    
    plt.suptitle(
        'Normalized sensitivity coefficients on integrated response (pc-Fos)',
        fontsize=24
    )
    plt.savefig('figure/sensitivity_h_PcFos_egf.pdf',bbox_inches='tight')
    plt.close()