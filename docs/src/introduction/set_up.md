## Prerequisites
In order to install and run BioMASS, a few preparations are required. To use BioMASS, you need to have Python3 installed, as well as a few packages.

The BioMASS repository is available on GitHub.

Install Git, open a terminal and run: 

`
$ sudo apt-get install git
`
<br>

Next, install Anaconda, a package manager, in order to handle dependencies such as matplotlib, numpy, etc.
Install anaconda3 and follow the prompts to proceed with the installation. <br>
Open a terminal and run:
`
$ wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
$ bash Anaconda3-2019.10-Linux-x86_64.sh
`

## Installation
The BioMASS repository can be installed directly from GitHub using the following command:

`
$ git clone https://github.com/okadalabipr/biomass.git
`
<br>

You can choose where to save the BioMASS repository.
After cloning, open your file browser and ensure that a copy of the BioMASS repository has been placed at your intended destination.

## Quick start

Follow these steps using command line:
If you quickly want to analyze a model which has already been implemented, open a terminal and navigate to the BioMASS repository.

We will use the "mapk_cascade" model for a quick parameter optimization and visualization of the results.

Start Python:

`$ python3
`

Load the model and the script required for optimization: <br>
`
$ >>> from biomass.models import mapk_cascade
`
`
$ >>> from biomass import optimize
`

Start the optimization by running the following commands:
`
$ >>> optimize(mapk_cascade, 1)
`
This step may take a few minutes. After finishing, one parameter set has been generated. Please check by navigating to the mapk_cascade directory. A new directory called "out" should have been created, containing a folder called "1". This folder contains the optimization results for the optimized parameter set 1. <br>

Let's visualize the optimized parameters. <br>
Type: <br>
`
$ >>>
`
<br>
We will explain the model files and functionality in more detail in the next steps.
