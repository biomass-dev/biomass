
# BioMASS tutorial

Biomass is a tool for parameter optimization. Users can supply BioMASS with differential equations, experimental conditions and experimental results to fit the equations to the experimental results. BioMASS can also perform sensitivity analysis to understand important components in complex biological network.

This tutorial is aimed for beginners and covers the basic setup required to run BioMASS as well as the steps needed to run BioMASS for custom differential equations.

This tutorial was made for the Laboratory of Cell Systems, Institute for Protein Research, Suita, Osaka, Japan.

For this tutorial, it will be helpful if users have basic background on bash command line.

## Machine setup
#
This step is needed to install the required softwares to run biomass

Install git for installing biomass.

```bash
$ sudo apt-get install git
```
<br>
Download biomass using git.

```bash
$ git clone https://github.com/okadalabipr/biomass.git
```

<br>
Install anaconda3 (Update based on the anaconda distribution).
Anaconda is a package manager needed for dependencies such as matplotlib, numpy, etc.
Follow the prompts to accept the installation.

```bash
$ wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
$ bash Anaconda3-2019.10-Linux-x86_64.sh
```
<br>
To test if python is installed, run the following line.

```bash
which python
```
<br>
Congrats! You should have a working python.

<br>


## BioMASS setup
#
There are several files that needs editing prior to running the program: 
<br>

1. Differential equations <br> 
```biomass/model/differential_equation.py```

2. Initial values <br>
```biomass/model/initial_condition.py```

3. Parameters <br>
```biomass/model/param_const.py```

4. Names of parameters and state variables (automated using bash script)<br>
```biomass/model/name2idx/parameters.py```<br>
```biomass/model/variables.py```

5. Functions to get the difference between experimental data and simulation<br>
```biomass/param_estim/fitness.py```

6. Observables, experimental data<br>
```biomass/biomass/observable.py```



<br>
To edit these files, go to the Biomass directory. We will be using the paper

[Oppelt *et al.*, ***NPJ Syst Biol***, 2018](https://doi.org/10.1038/s41540-018-0058-z), the equation is on the supplementary data p. 6.

<br>

### 1. Edit differential equations:
```biomass/model/differential_equation.py```

```python
from .name2idx import parameters as C
from .name2idx import variables as V

def diffeq(t,y,x):

    dydt = [0]*V.len_f_vars

    dydt[V.TNFR]     = 1*x[C.uptake]*x[C.TNF] -1*x[C.deact_TNFR]*y[V.TNFR]
    dydt[V.Ikk]      = -1*(x[C.act_Ikk_by_TNF]*y[V.TNFR]*y[V.Ikk]) +1*(x[C.trigger_iIkk]*y[V.iIkk])
    dydt[V.pIkk]     = 1*(x[C.act_Ikk_by_TNF]*y[V.TNFR]*y[V.Ikk]) -1*(x[C.act_pIkk]*y[V.pIkk])
    dydt[V.ppIkk]    = 1*(x[C.act_pIkk]*y[V.pIkk]) -1*(x[C.deact_ppIkk]*y[V.ppIkk])
    dydt[V.iIkk]     = 1*(x[C.deact_ppIkk]*y[V.ppIkk]) -1*(x[C.trigger_iIkk]*y[V.iIkk])
    dydt[V.NfkIkb]   = -1*((x[C.act_Ikb_by_Ikk])*y[V.pIkk]*y[V.NfkIkb]) -1*((x[C.act_Nfk_by_Ikk])*y[V.pIkk]*
                       y[V.NfkIkb]) +1*(x[C.form_complex]*y[V.Nfk]*y[V.Ikb]) +1*(x[C.ext_nNfkIkb]*y[V.nNfkIkb])\
                       *(x[C.Vnuc]/1)
    dydt[V.NfkpIkb]  = 1*((x[C.act_Ikb_by_Ikk])*y[V.pIkk]*y[V.NfkIkb]) -1*((x[C.act_Nfk_by_Ikk_complex])*
                       y[V.pIkk]*y[V.NfkpIkb]) -1*(x[C.split_NfkpIkb]*y[V.NfkpIkb])
    dydt[V.pNfkIkb]  = -1*((x[C.act_Ikb_complex])*y[V.pNfkIkb]) -1*((x[C.act_Ikb_by_Ikk])*y[V.pIkk]*
                       y[V.pNfkIkb]) +1*((x[C.act_Nfk_by_Ikk])*y[V.pIkk]*y[V.NfkIkb])
    dydt[V.pNfkpIkb] = 1*((x[C.act_Ikb_complex])*y[V.pNfkIkb]) +1*((x[C.act_Ikb_by_Ikk])*y[V.pIkk]*
                       y[V.pNfkIkb]) +1*((x[C.act_Nfk_by_Ikk_complex])*y[V.pIkk]*y[V.NfkpIkb]) -1*\
                       (x[C.split_NfkIkb]*y[V.pNfkpIkb])
    dydt[V.pNfk]     = 1*(x[C.split_NfkIkb]*y[V.pNfkpIkb]) -1*((x[C.int_Nfk]*x[C.eta_int_pNfk])*y[V.pNfk])
    dydt[V.Nfk]      = 1*(x[C.split_NfkpIkb]*y[V.NfkpIkb]) -1*(x[C.form_complex]*y[V.Nfk]*y[V.Ikb]) -1*\
                       ((x[C.int_Nfk])*y[V.Nfk])
    dydt[V.pIkb]     = 1*(x[C.split_NfkpIkb]*y[V.NfkpIkb]) +1*(x[C.split_NfkIkb]*y[V.pNfkpIkb]) -1*\
                       (x[C.degrad_Ikb]*y[V.pIkb])
    dydt[V.Ikb]      = -1*(x[C.form_complex]*y[V.Nfk]*y[V.Ikb]) +1*(x[C.prod_Ikb]*y[V.mIkb]) -1*(x[C.int_Ikb]*
                       y[V.Ikb])
    dydt[V.mIkb]     = 1*(x[C.prod_mIkb_by_nNfk]*y[V.nNfk]) -1*(x[C.degrad_mIkb]*y[V.mIkb])
    dydt[V.nIkb]     = 1*(x[C.int_Ikb]*y[V.Ikb])*(1/x[C.Vnuc]) -1*(x[C.form_complex_nuc]*y[V.nNfk]*y[V.nIkb])
    dydt[V.pnNfk]    = 1*((x[C.int_Nfk]*x[C.eta_int_pNfk])*y[V.pNfk])*(1/x[C.Vnuc]) -1*(x[C.deact_pnNfk]*
                       y[V.pnNfk])
    dydt[V.nNfk]     = 1*((x[C.int_Nfk])*y[V.Nfk])*(1/x[C.Vnuc]) +1*(x[C.deact_pnNfk]*y[V.pnNfk]) -1*\
                       (x[C.form_complex_nuc]*y[V.nNfk]*y[V.nIkb])
    dydt[V.nNfkIkb]  = 1*(x[C.form_complex_nuc]*y[V.nNfk]*y[V.nIkb]) -1*(x[C.ext_nNfkIkb]*y[V.nNfkIkb])
    dydt[V.RnaA20_1] = 1*(x[C.build_RnaA20]*y[V.nNfk]) -1*(x[C.shuttle_RnaA20]*y[V.RnaA20_1])
    dydt[V.RnaA20]   = 1*(x[C.shuttle_RnaA20]*y[V.RnaA20_1]) -1*(x[C.degrad_RnaA20]*y[V.RnaA20])
    dydt[V.A20]      = 1*(x[C.build_A20]*y[V.RnaA20]) -1*(x[C.degrad_A20]*y[V.A20])

    return dydt

```
<br>

### 2. Edit initial values: 
```biomass/model/initial_condition.py```

```python
from .name2idx import variables as V

def initial_values():
    y0 = [0]*V.len_f_vars

    y0[V.TNFR]     = 0.
    y0[V.Ikk]      = 1.
    y0[V.pIkk]     = 0.
    y0[V.ppIkk]    = 0.
    y0[V.iIkk]     = 0.
    y0[V.NfkIkb]   = 1.
    y0[V.NfkpIkb]  = 0.
    y0[V.pNfkIkb]  = 0.
    y0[V.pNfk]     = 0.
    y0[V.Nfk]      = 0.
    y0[V.pIkb]     = 0.
    y0[V.Ikb]      = 0.
    y0[V.mIkb]     = 0.
    y0[V.nIkb]     = 0.
    y0[V.pnNfk]    = 0.
    y0[V.nNfk]     = 0.
    y0[V.nNfkIkb]  = 0.
    y0[V.RnaA20_1] = 0.
    y0[V.RnaA20]   = 0.
    y0[V.A20]      = 0.

    return y0
```
<br>

### 3. Edit model parameters (From table S9): 
```biomass/model/param_const.py```

Make sure to change the line ```def f_params_noDCF():``` to ```def f_params:``` 

```python
from .name2idx import parameters as C

def f_params():
    x = [0]*C.len_f_params

    x[C.uptake] = 1.0000
    x[C.TNF] = 1.0000
    x[C.trigger_iIkk] = 0.0041
    x[C.deact_TNFR] = 0.0010
    x[C.deact_ppIkk] = 0.1660
    x[C.deact_pnNfk] = 1000.0000
    x[C.act_Ikk_by_TNF] = 0.0714
    x[C.act_pIkk] = 0.0648
    x[C.act_Ikb_by_Ikk] = 0.3980
    x[C.act_Nfk_by_Ikk] = 0.6438
    x[C.act_Nfk_by_Ikk_complex] = 0.2816
    x[C.act_Ikb_complex] = 1.3897
    x[C.form_complex] = 2.8390
    x[C.form_complex_nuc] = 1000.0000
    x[C.ext_nNfkIkb] = 1000.0000
    x[C.Vnuc] = 1.0000
    x[C.split_NfkpIkb] = 0.0811
    x[C.split_NfkIkb] = 1.0000
    x[C.int_Nfk] = 0.0100
    x[C.int_Ikb] = 0.1226
    x[C.eta_int_pNfk] = 17.9585
    x[C.degrad_Ikb] = 0.6308
    x[C.degrad_mIkb] = 0.0313
    x[C.degrad_RnaA20] = 0.0089
    x[C.degrad_A20] = 0.0116
    x[C.prod_Ikb] = 1.0000
    x[C.prod_mIkb_by_nNfk] = 0.0047
    x[C.build_RnaA20] = 1.0000
    x[C.build_A20] = 0.0006
    x[C.shuttle_RnaA20] = 0.0311

    return x
```
<br>

### 4. Use bash script to produce file of variables and parameters

Then, run the bash script to produce file of variables and parameters used in the analysis.

```biomass/model/name2idx/parameters.py```
<br>
```biomass/model/name2idx/variables.py```


```bash
# navigate to the tutorial folder 
$ cd biomass/tutorial
# run script
$ sh mk_param_var.sh
```
<br>

After producing the variables needed. We also need to specify the parameters for optimization. <br>
```biomass/param_estim/search_parameter.py``` <br>
This can also be automated by using the following bash script:

```bash
# navigate to the tutorial folder 
$ cd biomass/tutorial
# run script
$ sh mk_param_search.sh
```
<br>

After running the script, you will have a new ```search_parameter.py``` which contains the parameters you have specified in ```param_const.py```.

If you only need to optimize several parameters, you can open the ```search_parameter.py``` and **remove** the unwanted lines containing the parameters you don't want. The following lines will represent the parameter indices for optimization.

```python
		C.uptake,
		C.TNF,
		C.trigger_iIkk,
		C.deact_TNFR,
		C.deact_ppIkk,
		C.deact_pnNfk,
		C.act_Ikk_by_TNF,
		C.act_pIkk,
		C.act_Ikb_by_Ikk,
		C.act_Nfk_by_Ikk,
		C.act_Nfk_by_Ikk_complex,
		C.act_Ikb_complex,
		C.form_complex,
		C.form_complex_nuc,
		C.ext_nNfkIkb,
		C.Vnuc,
		C.split_NfkpIkb,
		C.split_NfkIkb,
		C.int_Nfk,
		C.int_Ikb,
		C.eta_int_pNfk,
		C.degrad_Ikb,
		C.degrad_mIkb,
		C.degrad_RnaA20,
		C.degrad_A20,
		C.prod_Ikb,
		C.prod_mIkb_by_nNfk,
		C.build_RnaA20,
		C.build_A20,
		C.shuttle_RnaA20,
```
<br>

Furthermore, if you want to specify the search region for different parameters, you can scroll down and uncomment the following lines:

```python
    # search_region[:, C.n10] = [1.00, 4.00]
    # search_region[:, C.n31] = [1.00, 4.00]
    # search_region[:, C.n57] = [1.00, 4.00]
    # search_region[:, C.nF31] = [1.00, 4.00]
```
<br>

You can then edit the parameter name and set the lower and upper limits.

```python
    search_region[:, C.uptake] = [0.5, 1.5]
    # search_region[:, C.n31] = [1.00, 4.00]
    # search_region[:, C.n57] = [1.00, 4.00]
    # search_region[:, C.nF31] = [1.00, 4.00]
```
<br>

### 5. Edit ```fitness.py``` file

When the simulations don't seem to work well, it might be better to change the function used to calculate the difference between the simulation and the experimental data. The default used is RSS and you can change it to cosine distance by changing the following line (line 61) from:

```python
                error[i] = compute_objval_rss(
```
<br>

to:

```python
                error[i] = compute_objval_cos(
```
<br>


Furthermore, when there is no constraints assuming that some parameters are equal to some other parameters, the lines containing constraints can be removed.

Remove/uncomment the following lines

from:
```python
    # constraints --------------------------------------------------------------
    x[C.V6] = x[C.V5]
    x[C.Km6] = x[C.Km5]
    x[C.KimpDUSP] = x[C.KimDUSP]
    x[C.KexpDUSP] = x[C.KexDUSP]
    x[C.KimpcFOS] = x[C.KimFOS]
    x[C.KexpcFOS] = x[C.KexFOS]
    x[C.p52] = x[C.p47]
    x[C.m52] = x[C.m47]
    x[C.p53] = x[C.p48]
    x[C.p54] = x[C.p49]
    x[C.m54] = x[C.m49]
    x[C.p55] = x[C.p50]
    x[C.p56] = x[C.p51]
    x[C.m56] = x[C.m51]
    # --------------------------------------------------------------------------
```

to:
```python
    # constraints --------------------------------------------------------------
    '''x[C.V6] = x[C.V5]
    x[C.Km6] = x[C.Km5]
    x[C.KimpDUSP] = x[C.KimDUSP]
    x[C.KexpDUSP] = x[C.KexDUSP]
    x[C.KimpcFOS] = x[C.KimFOS]
    x[C.KexpcFOS] = x[C.KexFOS]
    x[C.p52] = x[C.p47]
    x[C.m52] = x[C.m47]
    x[C.p53] = x[C.p48]
    x[C.p54] = x[C.p49]
    x[C.m54] = x[C.m49]
    x[C.p55] = x[C.p50]
    x[C.p56] = x[C.p51]
    x[C.m56] = x[C.m51]'''
    # --------------------------------------------------------------------------
```
<br>

### 6. Edit observables file<br>
Finally, you have to input the properties of the simulation, experimental information and data onto the observables file. <br>
```biomass/biomass/observable.py```

There are multiple lines that needs editing and this varies depending on the model.

* First change the observables. As the observable in this case is only the nuclear NF-kB, change the lines to the following lines from:

```python
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
```

to:
```python
observables = [
    'Nuclear_NFkB',
]
```
<br>

* We also need to change the time span, which is the total time span of the experiments. In this case, the data is from 0 to 5400 seconds, so, change the following from:
```python
    tspan = [0, 180]  # [start, end] (Unit time: 1 sec.)
```

to:
```python
    tspan = [0, 5400]  # [start, end] (Unit time: 1 sec.)
```
<br>

* To change the experimental conditions, edit the following lines based on the condition that the stimulus is only TNF from:

```python
    # Experimental conditions
    conditions = ['EGF', 'HRG']
```
to:
```python
    # Experimental conditions
    conditions = ['TNF']
```
<br>

* We need to change the simulation, first remove/comment out the steady state. In this example, the steady state was not simulated before TNF stimuli. 

from:
```python
    def simulate(self, x, y0):
        # get steady state
        x[C.Ligand] = x[C.no_ligand]  # No ligand
        (T_steady_state, Y_steady_state) = self._get_steady_state(
            diffeq, y0, self.tspan, tuple(x)
        )
        if T_steady_state < self.tspan[-1]:
            return False
        else:
            y0 = Y_steady_state[:]
```


to:
```python
    def simulate(self, x, y0):
        '''
        # get steady state
        x[C.TNF] = 0.0  # No ligand
        (T_steady_state, Y_steady_state) = self._get_steady_state(
            diffeq, y0, self.tspan, tuple(x)
        )
        if T_steady_state < self.tspan[-1]:
            return False
        else:
            y0 = Y_steady_state[:]
        '''
```
<br>

* Then, change the ligand input from EGF and HRG to TNF only, where the TNF is equal to the value of TNF specified in ```param_const.py```

from:
```python
        # add ligand
        for i, condition in enumerate(self.conditions):
            if condition == 'EGF':
                x[C.Ligand] = x[C.EGF]
            elif condition == 'HRG':
                x[C.Ligand] = x[C.HRG]
```

to:
```python
        # add ligand
        for i, condition in enumerate(self.conditions):
            if condition == 'TNF':
                x[C.TNF] = 1.0
```


* Change the observables index respectively. 
The following equation means that NF-kB observable includes the total of phosphorylated nuclear NF-kB ```V.pnNfk```, nuclear NF-kB ```V.nNfk```, and NF-kB IkB complex ```V.nNfkIkb```  . <br>
``` Y[:, V.pnNfk] + Y[:, V.nNfk] + Y[:, V.nNfkIkb]```

from:
```python
            (T, Y) = self._solveode(diffeq, y0, self.tspan, tuple(x))

            if T[-1] < self.tspan[-1]:
                return False
            else:
                self.simulations[observables.index('Phosphorylated_MEKc'), :, i] = (
                    Y[:, V.ppMEKc]
                )
                self.simulations[observables.index('Phosphorylated_ERKc'), :, i] = (
                    Y[:, V.pERKc] + Y[:, V.ppERKc]
                )
                self.simulations[observables.index('Phosphorylated_RSKw'), :, i] = (
                    Y[:, V.pRSKc] + Y[:, V.pRSKn]*(x[C.Vn]/x[C.Vc])
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
                    (Y[:, V.pcFOSn] + Y[:, V.cFOSn])*(x[C.Vn]/x[C.Vc])
                    + Y[:, V.cFOSc] + Y[:, V.pcFOSc]
                )
                self.simulations[observables.index('Phosphorylated_cFos'), :, i] = (
                    Y[:, V.pcFOSn]*(x[C.Vn]/x[C.Vc]) + Y[:, V.pcFOSc]
                )
```

to:
```python
            (T, Y) = self._solveode(diffeq, y0, self.tspan, tuple(x))

            if T[-1] < self.tspan[-1]:
                return False
            else:
                self.simulations[observables.index('Nuclear_NFkB'), :, i] = (
                    Y[:, V.pnNfk] + Y[:, V.nNfk] + Y[:, V.nNfkIkb]
                )
                
```
<br>

* Edit the time points in the experiment. As there is only one observable, only one set of time points will be required.

from:
```python
class ExperimentalData(object):

    experiments = [None]*len(observables)

    t2 = [0, 300, 600, 900, 1800, 2700, 3600, 5400]

    experiments[observables.index('Phosphorylated_MEKc')] = {
        'EGF': [0.000, 0.773, 0.439, 0.252, 0.130, 0.087, 0.080, 0.066],
        'HRG': [0.000, 0.865, 1.000, 0.837, 0.884, 0.920, 0.875, 0.789],
    }
    experiments[observables.index('Phosphorylated_ERKc')] = {
        'EGF': [0.000, 0.867, 0.799, 0.494, 0.313, 0.266, 0.200, 0.194],
        'HRG': [0.000, 0.848, 1.000, 0.971, 0.950, 0.812, 0.747, 0.595],
    }
    experiments[observables.index('Phosphorylated_RSKw')] = {
        'EGF': [0, 0.814, 0.812, 0.450, 0.151, 0.059, 0.038, 0.030],
        'HRG': [0, 0.953, 1.000, 0.844, 0.935, 0.868, 0.779, 0.558],
    }
    experiments[observables.index('Phosphorylated_cFos')] = {
        'EGF': [0, 0.060, 0.109, 0.083, 0.068, 0.049, 0.027, 0.017],
        'HRG': [0, 0.145, 0.177, 0.158, 0.598, 1.000, 0.852, 0.431],
    }

    # --------------------------------------------------------------------------
    t3 = [0, 600, 1800, 3600, 5400]

    experiments[observables.index('Phosphorylated_CREBw')] = {
        'EGF': [0, 0.446, 0.030, 0.000, 0.000],
        'HRG': [0, 1.000, 0.668, 0.460, 0.340],
    }

    # --------------------------------------------------------------------------
    t4 = [0, 600, 1200, 1800, 2700, 3600, 5400]

    experiments[observables.index('cfos_mRNA')] = {
        'EGF': [0, 0.181, 0.476, 0.518, 0.174, 0.026, 0.000],
        'HRG': [0, 0.353, 0.861, 1.000, 0.637, 0.300, 0.059],
    }

    # --------------------------------------------------------------------------
    t5 = [0, 900, 1800, 2700, 3600, 5400]

    experiments[observables.index('cFos_Protein')] = {
        'EGF': [0, 0.078, 0.216, 0.240, 0.320, 0.235],
        'HRG': [0, 0.089, 0.552, 0.861, 1.000, 0.698],
    }
    experiments[observables.index('dusp_mRNA')] = {
        'EGF': [0.000, 0.177, 0.331, 0.214, 0.177, 0.231],
        'HRG': [0.000, 0.221, 0.750, 1.000, 0.960, 0.934],
    }
```
to:
```python
class ExperimentalData(object):

    experiments = [None]*len(observables)

    t = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]

    experiments[observables.index('Nuclear_NFkB')] = {
        'TNF': [0.00, 0.40, 1.00, 0.63, 0.28, 0.14, 0.16, 0.45, 0.64, 0.29, 0.18, 0.26, 0.40],
    }
```
<br>

* Finally, you just have to change the following lines. In this case, we only have one time point corresponding to the result Nuclear NF-kB. Take care that as there is only one experimental time point corresponding to ```self.t```. When there are more than one experimental data, time points should be written according to the different results such as in previous section ```t1```, ```t2```, ```t3```, and so forth.

from:
```python
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
```
<br>

to:
```python
    def get_timepoint(self, obs_idx):
        exp_t = self.t

        return list(map(int, exp_t))
```



<br>
Now, the setup is finished and BioMASS can be run.

<br>
<br>

## Running BioMASS
#

### ***Running the simulations***

Firstly, move back to the home directory of BioMASS.
```bash
$ cd ..
```
<br>

Then we can run the following to run simulations. Please take note of the following points:
* ```nohup``` is not necessary but helpful so that the process can run when logged out. 
* To run a program in the background, enter the command for that job, followed by the ```&``` sign.
* ```python optimize.py``` is the optimization script
* ```n``` is the set of random  

For example, to run simulations with random parameter sets from ```1``` **to** ```5``` (```1```, ```2```, ```3```, ```4```, ```5```) you would want to write:
```bash
$ nohup python optimize.py 1 5 &
```
<br>

For only one random parameter set number ```4```.
```bash
$ nohup python optimize.py 4 &
```
<br>

### ***Visualization of simulation results***

Before running visualization, we have to firstly set the plot by editing ```/biomass/param_estim/plot_func.py```

Keep in mind that this current model has the time units in minutes and therefore the following lines have to be edited (line 78) (for models that uses seconds)

from:
```python
                    plt.plot(
                        np.array(exp_t) / 60., exp.experiments[i][condition],
                        shape[l], markerfacecolor='None', markeredgecolor=cmap[l],
                        clip_on=False
                    )
```

to:
```python
                    plt.plot(
                        np.array(exp_t), exp.experiments[i][condition],
                        shape[l], markerfacecolor='None', markeredgecolor=cmap[l],
                        clip_on=False
                    )
```
<br>

Also the following lines (x-axis parameters) have to be changed to accomodate the results of this model. (line 83-84)

from:
```python
        plt.xlim(0, 90)
        plt.xticks([0, 30, 60, 90])
```
<br>

to:
```python
        plt.xlim(0, 180)
        plt.xticks([0, 60, 120, 180])
```
<br>

To visualize the simulation results you can then use the function ```run_sim.py``` with several visualization options:

* Visualization type: 
    * ```average``` : average of the results with parameter sets in ```out/```
    * ```best``` : the best results (what is best?)
    * ```original```: original result before optimization
    * ```n```: simulation with the parameter set in ```out/n```

* Show all: ```show_all``` 

* Standard deviation: ```stdev``` to plot the standard deviation. (only available for ```average``` visualization type)

```bash
# to get the average, show all and visualize standard deviation by error bars
$ python run_sim.py best show_all stdev
```
<br>

### ***Sensitivity analysis***

You can calculate sensitivity coefficients on rate equations and non-zero initial values using ```analyze.py```

To obtain values for sensitivity of the rate equations, the time derivatives of state variables must be described via rate equations (See ```differential_equation.py```) and you need to edit ```biomass/analysis/reaction/reaction.py```

In this case, we just have to change the following (line 24) of ```analyze.py``` from:

```python
    reaction.sensitivity_barplot(metric=str(args[1]))
#   reaction.sensitivity_heatmap(metric='args[1]')
#   nonzero_init.sensitivity_barplot(metric=str(args[1]))
#   nonzero_init.sensitivity_heatmap(metric='args[1]')
```
<br>

to:
```python
#   reaction.sensitivity_barplot(metric=str(args[1]))
#   reaction.sensitivity_heatmap(metric='args[1]')
    nonzero_init.sensitivity_barplot(metric=str(args[1]))
#   nonzero_init.sensitivity_heatmap(metric='args[1]')
```
<br>

There are 3 options available:
* ```amplitude```: the maximum value.
* ```duration```: the time it takes to decline below 10% of its maximum.
* ```integral```: the integral of concentration over the observation time.

For example, to use the maximum value as a signaling metric:

```bash
$ python analyze.py amplitude
```
<br>

The final results should be available in the folder ```/figure```


If you have any questions, please email us at nico@protein.osaka-u.ac.jp or himoto@protein.osaka-u.ac.jp

## Author
- Johannes Nicolaus Wibisana  (2020/02/06)
