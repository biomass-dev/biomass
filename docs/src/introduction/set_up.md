## Prerequisites
In order to install and run BioMASS, a few preparations are required. To use BioMASS, you need an installation of Python3, as well as a few packages.

The BioMASS repository is available on GitHub.
If you have Python3 and Anaconda (or another package manager of your choice) already running, please proceed to the next section.

Linux:
To install Git, open a terminal and run: 

`
$ sudo apt-get install git
`
<br>

Next, install Anaconda, a package manager, in order to handle dependencies such as matplotlib, numpy, etc.
Install anaconda3 and follow the prompts to proceed with the installation. <br>
Open a terminal and run:
`
$ wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh 
`
<br>
`
$ bash Anaconda3-2019.10-Linux-x86_64.sh
`

## BioMASS Installation
The BioMASS repository can be installed directly from GitHub using the following command:

`
$ git clone https://github.com/okadalabipr/biomass.git
`
<br>

You can choose where to save the BioMASS repository.
After cloning, open your file browser and ensure that a copy of the BioMASS repository has been placed at your intended destination.

## Quick start

Follow these steps using command line:
If you quickly want to simulate and visualize a model which has already been implemented, open a terminal and navigate to the BioMASS repository.

We will use the "mapk_cascade" model for a quick parameter optimization and visualization of the results.

Start Python:

`$ python3
`

Load the model and the module required for simulation: <br>
`
$ >>> from biomass.models import mapk_cascade
`
<br>

`
$ >>> from biomass import run_simulation
`

Next, run the simulation: <br>
`
$ >>> run_simulation(mapk_cascade, viz_type='original', show_all=False, stdev=True)

`
The simulation uses the visualization type "original", meaning that it applies the parameters already defined in the model.
The simulation results will be saved in an automatically created directory "original" with the following destination:

<br>
`biomass->biomass-> models-> mapk_cascade -> figure -> simulation -> original`
<br>
Navigate to the `original` directory. Three files should have been created. Open them and check your simulation results! Solid lines correspond to simulated data, circles and triangles correspond to experimental data.

You can also test a parameter optimization. Please not that the following procedure may take up to one day to complete. 
Load the optimization module: <br>
`
$ >>> from biomass import optimize
`

Start the optimization for one parameter set: <br>
`
$ >>> optimize(mapk_cascade, 1)
`
After finishing, one parameter set has been generated. Please check by navigating to the mapk_cascade directory. A new directory called "out" should have been created, containing a folder called "1". This folder contains the optimization results for the optimized parameter set 1. <br>

Let's visualize the optimized parameters. <br>
Type: <br>
`
$ >>> from biomass import run_simulation
`
<br>
`
$ >>> from biomass.models import mapk_cascade 
`
<br>

`
$ >>> run_simulation(mapk_cascade, viz_type='best',  show_all=False, stdev=True)
`
<br>
Note: Both for parameter optimization and model simulation, your working directory needs to be the 'biomass' folder. In this example, we have only optimized one parameter set, and hence, have to choose the visualization type 'best'. You can change the visualization type to 'average' if you optimize more than one parameter set.
We will explain the model files and functionality in more detail in the next steps.
