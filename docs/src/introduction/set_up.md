## Prerequisites
In order to install and run BioMASS, a few preparations are required.


## Installation
The BioMASS repository can be installed directly from GitHub using the following command:

`
$ git clone https://github.com/okadalabipr/biomass.git
`
You can choose where to save the BioMASS repository.
After cloning, open your file browser and ensure that a copy of the BioMASS repository has been placed at your intended destination.

## Quick start

Follow these steps using command line:
If you quickly want to analyze a model which has already been implemented, open a terminal/console and navigate to the BioMASS repository.

We will use the "mapk_cascade" model for a quick parameter optimization and visualization of the results.

Start Python:

`$ python3
`

Load the model and the script required for optimization:
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

Let's visualize the optimized parameters.
Type:
`
$ >>>
`

We will explain the model files and functionality in more detail in the next steps.
