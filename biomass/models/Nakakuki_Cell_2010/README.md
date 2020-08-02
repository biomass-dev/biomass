# Nakakuki_Cell_2010
This repository contains data and modeling code for the following paper:

Nakakuki, T. *et al.* Ligand-specific c-Fos expression emerges from the spatiotemporal control of ErbB network dynamics. *Cell* **141**, 884–896 (2010). https://doi.org/10.1016/j.cell.2010.03.054

## BioModels
- [Nakakuki2010_CellFateDecision_Mechanistic](https://www.ebi.ac.uk/biomodels/BIOMD0000000250)

## Description
A brief description of each file is below:

|Name|Contents|
|---|---|
|[`name2idx/`](./name2idx/)|Names of model parameters and species|
|[`set_model.py`](./set_model.py)|Differential equation, parameters and initial condition|
|[`observalbe.py`](./observable.py)|Observables, simulations and experimental data|
|[`viz.py`](./viz.py)|Plotting parameters for customizing figure properties|
|[`set_search_param.py`](./set_search_param.py)|Model parameters to optimize and search region|
|[`fitness.py`](./fitness.py)|An objective function to be minimized, i.e., the distance between model simulation and experimental data|
|[`reaction_network.py`](./reaction_network.py)|Reaction indices grouped according to biological processes|