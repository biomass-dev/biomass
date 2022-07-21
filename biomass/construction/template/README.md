# Template for BioMASS model construction

A brief description of each file/folder is below:

| Name                                           | Content                                                                                                  |
| ---------------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| [`name2idx/`](./name2idx/)                     | Names of model parameters and species                                                                    |
| [`reaction_network.py`](./reaction_network.py) | Flux vector and reaction indices grouped according to biological processes                               |
| [`ode.py`](./ode.py)                           | Differential equation, parameters and initial condition                                                  |
| [`observalbe.py`](./observable.py)             | Observables, simulations and experimental data                                                           |
| [`viz.py`](./viz.py)                           | Plotting parameters for customizing figure properties                                                    |
| [`search_param.py`](./search_param.py)         | Lower and upper bounds of model parameters to be estimated                                               |
| [`problem.py`](./problem.py)                   | An objective function to be minimized, i.e., the distance between model simulation and experimental data |

Example models can be found at https://biomass-core.readthedocs.io/en/latest/models.html.