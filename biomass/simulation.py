import numpy as np
from scipy.integrate import ode

from .model.name2idx import parameters as C
from .model.name2idx import variables as V
from .model.differential_equation import diffeq
from .observable import num_observables,species

def solveode(diffeq,y0,tspan,args):
    sol = ode(diffeq)
    sol.set_integrator('vode',method='bdf',min_step=1e-8,with_jacobian=True)
    sol.set_initial_value(y0,tspan[0])
    sol.set_f_params(args)

    T = [tspan[0]]
    Y = [y0]

    while sol.successful() and sol.t < tspan[-1]:
        sol.integrate(sol.t+1.)
        T.append(sol.t)
        Y.append(sol.y)

    return np.array(T),np.array(Y)


class NumericalSimulation(object):

    tspan = range(5401) # Unit time: 1 sec.
    condition = 2

    t = np.array(tspan)/60. # sec. -> min. (plot_func.py)
    
    simulations = np.empty((num_observables,len(tspan),condition))

    def simulate(self,x,y0):

        for i in range(self.condition):
            if i==0:
                x[C.Ligand] = x[C.EGF]
            elif i==1:
                x[C.Ligand] = x[C.HRG]

            (T,Y) = solveode(diffeq,y0,self.tspan,tuple(x))

            if T[-1] < self.tspan[-1]:
                return False
            else:
                self.simulations[species['Phosphorylated_MEKc'],:,i] = \
                    Y[:,V.ppMEKc]
                    
                self.simulations[species['Phosphorylated_ERKc'],:,i] = \
                    Y[:,V.pERKc] + Y[:,V.ppERKc]
                    
                self.simulations[species['Phosphorylated_RSKw'],:,i] = \
                    Y[:,V.pRSKc] + Y[:,V.pRSKn]*(x[C.Vn]/x[C.Vc])
                    
                self.simulations[species['Phosphorylated_CREBw'],:,i] = \
                    Y[:,V.pCREBn]*(x[C.Vn]/x[C.Vc])
                    
                self.simulations[species['dusp_mRNA'],:,i] = \
                    Y[:,V.duspmRNAc]
                    
                self.simulations[species['cfos_mRNA'],:,i] = \
                    Y[:,V.cfosmRNAc]
                    
                self.simulations[species['cFos_Protein'],:,i] = \
                    (Y[:,V.pcFOSn] + Y[:,V.cFOSn])*(x[C.Vn]/x[C.Vc]) + Y[:,V.cFOSc] + Y[:,V.pcFOSc]
                    
                self.simulations[species['Phosphorylated_cFos'],:,i] = \
                    Y[:,V.pcFOSn]*(x[C.Vn]/x[C.Vc]) + Y[:,V.pcFOSc]