# BioMASS overview
In this section, we will familiarize ourselves with the structure and functionality of BioMASS.
We will use a model available in the BioMASS repository.

Let's navigate to the models directory in the BioMASS repository. <br>
Go to: 
<br>
`
biomass -> biomass -> models 
`
<br>
Three models are available in the models directory, each stored in its own directory. 
Each of these model directories contains all files necessary to define a model in BioMASS. <br>

We will use the model of a MAPK cascade in this overview.

Open the mapk_cascade directory. <br>
`
models -> mapk_cascade
`
<br>
In the following, we will briefly introduce the relevant model files and their functions. In particular, we will focus on the following three files:
- `set_model.py` for defining the model itself
- `observable.py` for setting observables etc.
- `set_search_param.py` for defining the parameters to optimize


## *set_model.py*: Defining the model 

The file 'set_model.py' contains the structure of the model implemented as ODEs, as well as the definitions for all parameters and initial values. Please open the file 'set_model.py'. We will explain the three main functions in this file.

- *diffeq* 
Have a look at the *diffeq* function  starting in line 10. Here, the rate laws and the ODEs are defined. As you can see, this model contains in total eight ODEs defined from lines 25 to 32. 

- *param_values* 
Next, look at the *param_values* function starting from line 37. Here, all parameters that are used in the ODEs are defined.

- *initial_values*
Each ODE requires some initial value. These are defined in the function *initial_values* starting from line 66.
<br>
Note that you can choose which of the parameters and initial values will be optimized by editing the file 'set_search_param.py'.

## *observable.py*: Defining observables, adding experimental data, defining conditions, setting *in silico* simulation time
Open the file 'observable.py'. In this file, you can define which of the species contained in the model you would like to observe, that is, retrieving the dynamics of the simulated model, as well as add experimental data. Optionally, you can also define experimental conditions which you would like to explore "in silico." 

### Observables
As you can see, in this model, two observables are defined, the double phosphorylated MAPK (*biphosphorylated_MAPK*), and the unphosphorylated MAPK (*unphosphorylated_MAPK*). 
<br>
Starting from line 60, you need to specify which species in your model correspond to the observables. In this model, the two observables are defined by only one model species, each. However, you can add as many model species to your observable as required. For example, if you would like to observe the total amount of a molecule of interest, you can add the corresponding model species.
<br>
Note that for each observable defined from line 8, you need to specify the corresponding model species starting from line 58.

### Experimental data
The function *set_data* starting from line 172 contains the experimental data which is used by BioMASS for parameter optimization. In this model, we have experimental data points for the two observables which were defined earlier.

### Conditions
It is common that experiments are conducted in different conditions, such as the addition of a stimulant or inihibitor. You can recapulate such conditions with BioMASS. In this model, only a control condition is defined. 

### *in silico* simulation time
Lastly, in the function *get_timepoint* you can define the *in silico* time points for your simulation output. In this model, we observe from 0 to 150 with increments of 5 (line 202). Note that you can define the time period you want to study and the time resolution according to your own needs. The *in silico* does not have to match your experimental data points.

## *set_search_param.py*: Defining paramters to optimize


