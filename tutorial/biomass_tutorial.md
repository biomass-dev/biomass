
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
There are several files that needs editing prior to running the program: 

```bash
$ cd biomass/models/[your_model]
```
### 1. Edit names of model parameters and species:
- ```name2idx/```

    - ```parameters.py```
        - Write names of model parameters in ```parameters``` as string

    - ```species.py```
        - Write names of model species in ```species``` as string

- Or you can use bash script to produce file of parameters and species after editing ```set_model.py```.
    ```bash
    # navigate to the tutorial folder 
    $ cd biomass/tutorial
    # run script
    $ sh mk_param_var.sh
    ```

### 2. Edit differential equations, paramters, initial condition:
- ```set_model.py```

    - ```diffeq(t, y, x)```
        - Computes the derivative of ```y``` at ```t```.

    - ```param_values()```
        - Return model parameters.

    - ```initial_values()```
        - Initial conditions, specified as a vector. y0 must be the same length as the vector output of ```diffeq```, so that y0 contains an initial condition for each equation defined in ```diffeq```.

### 3. Edit observables, simulations and experimental data

- ```observable.py```

    This is the file to define the simulations you want to run and input the experimental data that you are going to use to try and fit the parameters to.


### 4. Edit an objective function for parameter estimation
- ```fitness.py```

    - ```objective(indiv_gene)``` : An objective function to be minimized, i.e., the distance between model simulation and experimental data 

    ```python
    def objective(indiv_gene):
        """Define an objective function to be minimized
        """
        indiv = decode_gene2val(indiv_gene)

        (x, y0) = update_param(indiv)

        exp = ExperimentalData()
        sim = NumericalSimulation()
        
        error = np.zeros(14)
        # sum of squared differences between the experimental data and the simulated values 

        norm_max = np.max(sim.simulations[observables.index('Phosphorylated_ERKc')])
        error[0] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_ERKc'), exp.t2, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_ERKc')]['EGF']
        )
        error[1] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_ERKc'), exp.t2, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_ERKc')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('Phosphorylated_RSKw')])
        error[2] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_RSKw'), exp.t2, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_RSKw')]['EGF']
        )
        error[3] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_RSKw'), exp.t2, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_RSKw')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('Phosphorylated_CREBw')])
        error[4] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_CREBw'), exp.t3, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_CREBw')]['EGF']
        )
        error[5] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_CREBw'), exp.t3, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_CREBw')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('dusp_mRNA')])
        error[6] = _compute_objval_rss(
            sim.simulations[observables.index('dusp_mRNA'), exp.t5, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('dusp_mRNA')]['EGF']
        )
        error[7] = _compute_objval_rss(
            sim.simulations[observables.index('dusp_mRNA'), exp.t5, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('dusp_mRNA')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('cfos_mRNA')])
        error[8] = _compute_objval_rss(
            sim.simulations[observables.index('cfos_mRNA'), exp.t4, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('cfos_mRNA')]['EGF']
        )
        error[9] = _compute_objval_rss(
            sim.simulations[observables.index('cfos_mRNA'), exp.t4, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('cfos_mRNA')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('cFos_Protein')])
        error[10] = _compute_objval_rss(
            sim.simulations[observables.index('cFos_Protein'), exp.t5, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('cFos_Protein')]['EGF']
        )
        error[11] = _compute_objval_rss(
            sim.simulations[observables.index('cFos_Protein'), exp.t5, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('cFos_Protein')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('Phosphorylated_cFos')])
        error[12] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_cFos'), exp.t2, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_cFos')]['EGF']
        )
        error[13] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_cFos'), exp.t2, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_cFos')]['HRG']
        )
    ```


### 5. Specify parameters to optimize and search region
- ```set_search_param.py```

    - ```get_search_index```

        - Specify names of model parameters and/or initial values to optimize.

    - ```get_search_region```

        - Set search region in optimization.

            - for all paramters

            ```python
            # search_idx[0] : model parameters
            for i, j in enumerate(search_idx[0]):
                search_rgn[0, j] = search_param[i] * 0.1  # lower bound
                search_rgn[1, j] = search_param[i] * 10.  # upper bound
            ```

            - for specific parameters

            ```python
            # Hill coefficient
            search_rgn[:, C.n10] = [1.00, 4.00]
            search_rgn[:, C.n31] = [1.00, 4.00]
            search_rgn[:, C.n57] = [1.00, 4.00]
            search_rgn[:, C.nF31] = [1.00, 4.00]
            ```

    - ```update_param```

        - To impose parameter value constraints,

        ```python
        def update_param(indiv):
            x = param_values()
            y0 = initial_values()

            search_idx = get_search_index()

            for i, j in enumerate(search_idx[0]):
                x[j] = indiv[i]
            for i, j in enumerate(search_idx[1]):
                y0[j] = indiv[i+len(search_idx[0])]

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

            return x, y0
        ```

### 6. Group reactions according to biological processes
- ```reaction_network.py```

    - For visualizing the result of sensitivity analysis

### 7. Make your model executable
- ```biomass/current_model.py```

    ```python
    # from biomass.models.Nakakuki_Cell_2010 import *
    from biomass.models.[your_model] import *
    ```

## Running BioMASS

### ***Parameter eatimation***

Firstly, move back to the home directory of BioMASS.
```bash
$ cd ..
```
<br>

Then we can run the following to run optimization. Please take note of the following points:
* ```nohup``` is not necessary but helpful so that the process can run when logged out. 
* To run a program in the background, enter the command for that job, followed by the ```&``` sign.
* ```python optimize.py``` is the optimization script.

For only one parameter set number ```1```,
```bash
$ nohup python optimize.py 1 &
```
<br>

For example, to run optimization with parameter sets from ```1``` **to** ```5``` (```1```, ```2```, ```3```, ```4```, ```5```) simultaneously,
```bash
$ nohup python optimize.py 1 5 &
```
<br>


### ***Visualization of simulation results***

To visualize the simulation results you can then use the function ```simulate_all``` with several visualization options:

* **Visualization type**: 
    * ```'average'``` : average of the results with parameter sets in ```out/```
    * ```'best'``` : the best results (what is best?)
    * ```'original'```: original result before optimization
    * ```'n'```: simulation with the parameter set in ```out/n```

* **Show all** : ```show_all``` 

* **Standard deviation** : ```stdev``` to plot the standard deviation. (only available for ```average``` visualization type)

```python
import warnings
warnings.filterwarnings('ignore')

from biomass.param_estim import simulate_all

# to get the average and visualize standard deviation by error bars
simulate_all(viz_type='average', show_all=False, stdev=True)
```
<br>

### ***Sensitivity analysis***

You can calculate sensitivity coefficients on rate equations and non-zero initial values using ```analyze```.

To obtain values for sensitivity of the rate equations, the time derivatives of state variables must be described via rate equations (See ```set_model.py```) and you need to edit ```reaction_network.py```

**metric**: 3 options available
- ```'amplitude'```
    : The maximum value.
- ```'duration'```
    : The time it takes to decline below 10% of its maximum.
- ```'integral'```
    : The integral of concentration over the observation time.

**style**: 2 options available
- ```'barplot'``` : To visualizehe the averaged sensitivity coefficients. 
- ```'heatmap'``` : To visualize the individual sensitivity coefficients for each parameter set.

For example, to calculate sensitivity coefficients on rate equations, use the maximum value as a signaling metric and save barplot:

```python
from biomass.analysis import reaction, nonzero_init

reaction.analyze(metric='amplitude', style='barplot')
```
<br>

The final results should be available in the folder ```/figure```


If you have any questions, please email us at nico@protein.osaka-u.ac.jp or himoto@protein.osaka-u.ac.jp

## Author
- Johannes Nicolaus Wibisana  (2020/02/06)
