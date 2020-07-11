# Nakakuki_Cell_2010
This repository contains data and modeling code for the following paper:

Nakakuki, T. *et al.* Ligand-specific c-Fos expression emerges from the spatiotemporal control of ErbB network dynamics. *Cell* **141**, 884–896 (2010). https://doi.org/10.1016/j.cell.2010.03.054

---
This mechanistic model describes the activation of immediate early genes such as cFos after epidermal growth factor (EGF) or heregulin (HRG) stimulation of the MAPK pathway. Phosphorylated cFos is a key transcription factor triggering downstream cascades of cell fate determination. The model can explain how the switch-like response of p-cFos emerges from the spatiotemporal dynamics. This mechanistic model comprises the explicit reaction kinetics of the signal transduction pathway, the transcriptional and the posttranslational feedback and feedforward loops. In the article, two different mechanistic models have been studied, the first one based on previously known interactions but failing to account for the experimental data and the second one including additional interactions which were discovered and confirmed by new experiments. The mechanistic model encoded here is the second one, the extended and at the time of creation most complete model of cell fate decision making in response to different doses of EGF or HRG stimulation.

## Description
A brief description of each file is below:

|Name|Contents|
|---|---|
|[`name2idx/`](./name2idx/)|Names of model parameters and species|
|[`set_model.py`](./set_model.py)|Differential equation, parameters and initial condition|
|[`observalbe.py`](./observable.py)|observables, simulations, experimental data and plotting options|
|[`set_search_param.py`](./set_search_param.py)|Model parameters to optimize and search region|
|[`fitness.py`](./fitness.py)|An objective function to be minimized, i.e., the distance between model simulation and experimental data|
|[`reaction_network.py`](./reaction_network.py)|Reaction indices grouped according to biological processes|