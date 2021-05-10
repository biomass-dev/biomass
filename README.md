# BioMASS

[![Actions Status](https://github.com/okadalabipr/biomass/workflows/Tests/badge.svg)](https://github.com/okadalabipr/biomass/actions)
[![PyPI version](https://img.shields.io/pypi/v/biomass.svg?logo=PyPI&logoColor=white)](https://pypi.python.org/pypi/biomass/)
[![License](https://img.shields.io/badge/License-Apache%202.0-green.svg)](https://opensource.org/licenses/Apache-2.0)
[![Downloads](https://pepy.tech/badge/biomass)](https://pepy.tech/project/biomass)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/biomass.svg)](https://pypi.python.org/pypi/biomass/)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/okadalabipr/biomass.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/okadalabipr/biomass/context:python)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

<img align="left" src="https://raw.githubusercontent.com/okadalabipr/biomass/master/resources/images/logo.png" width="300">

Mathematical modeling is a powerful method for the analysis of complex biological systems. Although there are many researches devoted on producing models to describe dynamical cellular signaling systems, most of these models are limited and do not cover multiple pathways. Therefore, there is a challenge to combine these models to enable understanding at a larger scale. Nevertheless, larger network means that it gets more difficult to estimate parameters to reproduce dynamic experimental data needed for deeper understanding of a system.

To overcome this problem, we developed _BioMASS_, a Python framework for **M**odeling and **A**nalysis of **S**ignaling **S**ystems. The BioMASS framework allows efficient optimization of multiple parameter sets simultaneously and generates the multiple parameter candidates that explain the signaling dynamics of interest. These parameter candidates can be further evaluated by their distribution and sensitivity analysis as a part of alternative information about the hidden regulatory mechanism of the system.

## Features

- Parameter estimation of ODE models
- Local sensitivity analysis
- Effective visualization of simulation results

## Installation

The BioMASS library is available on [PyPI](https://pypi.org/project/biomass/).

```bash
$ pip install biomass
```

BioMASS supports Python 3.7 or newer.

## Usage

- [Circadian clock](https://github.com/okadalabipr/biomass/tree/master/biomass/models/circadian_clock)
- [Immediate-early gene response](https://github.com/okadalabipr/biomass/tree/master/biomass/models/Nakakuki_Cell_2010)
- [Insulin signaling](https://github.com/okadalabipr/biomass/tree/master/biomass/models/insulin_signaling)
- [MAPK cascade](https://github.com/okadalabipr/biomass/tree/master/biomass/models/mapk_cascade)
- [NF-κB pathway](https://github.com/okadalabipr/biomass/tree/master/biomass/models/nfkb_pathway)
- [TGF-β/SMAD pathway](https://github.com/okadalabipr/biomass/tree/master/biomass/models/tgfb_smad)

We will use the model of immediate-early gene response ([Nakakuki_Cell_2010](https://github.com/okadalabipr/biomass/tree/master/biomass/models/Nakakuki_Cell_2010)) for parameter estimation, visualization of simulation results and sensitivity analysis.

### Model Preparation

A brief description of each file/folder is below:

| Name                  | Content                                                                                                  |
| --------------------- | -------------------------------------------------------------------------------------------------------- |
| `name2idx/`           | Names of model parameters and species                                                                    |
| `set_model.py`        | Differential equation, parameters and initial condition                                                  |
| `observalbe.py`       | Observables, simulations and experimental data                                                           |
| `viz.py`              | Plotting parameters for customizing figure properties                                                    |
| `set_search_param.py` | Model parameters to optimize and search region                                                           |
| `fitness.py`          | An objective function to be minimized, i.e., the distance between model simulation and experimental data |
| `reaction_network.py` | Reaction indices grouped according to biological processes                                               |

```python
from biomass import Model
from biomass.models import Nakakuki_Cell_2010

model = Model(Nakakuki_Cell_2010.__package__).create(show_info=True)
```

↓

```
Nakakuki_Cell_2010 information
------------------------------
36 species
115 parameters, of which 75 to be estimated
```

### Parameter Estimation of ODE Models (_n_ = 1, 2, 3, · · ·)

Parameters are adjusted to minimize the distance between model simulation and experimental data.

```python
from biomass import optimize

optimize(
    model=model, start=1, options={
        "popsize": 3,
        "max_generation": 100,
        "allowable_error": 0.5,
        "local_search_method": "DE",
        "maxiter": 50,
    }
)
```

The temporary result will be saved in **out/n/** after each iteration.

Progress list: **out/n/optimization.log**

```
Generation1: Best Fitness = 5.864228e+00
Generation2: Best Fitness = 5.864228e+00
Generation3: Best Fitness = 4.488934e+00
Generation4: Best Fitness = 3.793744e+00
Generation5: Best Fitness = 3.652047e+00
Generation6: Best Fitness = 3.652047e+00
Generation7: Best Fitness = 3.652047e+00
Generation8: Best Fitness = 3.452999e+00
Generation9: Best Fitness = 3.180878e+00
Generation10: Best Fitness = 1.392501e+00
Generation11: Best Fitness = 1.392501e+00
Generation12: Best Fitness = 1.392501e+00
Generation13: Best Fitness = 1.392501e+00
Generation14: Best Fitness = 7.018051e-01
Generation15: Best Fitness = 7.018051e-01
Generation16: Best Fitness = 7.018051e-01
Generation17: Best Fitness = 7.018051e-01
Generation18: Best Fitness = 7.018051e-01
Generation19: Best Fitness = 6.862063e-01
Generation20: Best Fitness = 6.862063e-01
```

- If you want to continue from where you stopped in the last parameter search,

```python
from biomass import optimize_continue

optimize_continue(
    model=model, start=1, options={
        "popsize": 3,
        "max_generation": 200,
        "allowable_error": 0.5,
        "local_search_method": "DE",
        "maxiter": 50,
    }
)
```

- If you want to search multiple parameter sets (e.g., from 1 to 10) simultaneously,

```python
from biomass import optimize

optimize(
    model=model, start=1, end=10, options={
        "popsize": 5,
        "max_generation": 100,
        "allowable_error": 0.5,
        "local_search_method": "DE",
        "maxiter": 50,
    }
)
```

- Data Export and Visualization

```python
from biomass.result import OptimizationResults

res = OptimizationResults(model)
# Export optimized parameters in CSV format
res.to_csv()
# Visualize objective function traces for different optimization runs
res.trace_obj()
```

### Visualization of Simulation Results

```python
from biomass import run_simulation

run_simulation(model, viz_type='average', show_all=False, stdev=True)
```

![simulation_average](https://raw.githubusercontent.com/okadalabipr/biomass/master/resources/images/simulation_average.png)

Points (blue diamonds, EGF; red squares, HRG) denote experimental data, solid lines denote simulations

### Sensitivity Analysis

The single parameter sensitivity of each reaction is defined by<br>

_s<sub>i</sub>_(_q_(**v**),_v<sub>i</sub>_) = _∂_ ln(_q_(**v**)) / _∂_ ln(_v<sub>i</sub>_) = _∂_ _q_(**v**) / _∂_ _v<sub>i</sub>_ · _v<sub>i</sub>_ / _q_(**v**)

where _v<sub>i</sub>_ is the _i_<sup>th</sup> reaction rate, **v** is reaction vector **v** = (_v<sub>1</sub>_, _v<sub>2</sub>_, ...) and _q_(**v**) is a target function, e.g., time-integrated response, duration. Sensitivity coefficients were calculated using finite difference approximations with 1% changes in the reaction rates.

```python
from biomass import run_analysis

run_analysis(model, target='reaction', metric='integral', style='barplot')
```

![sensitivity_PcFos](https://raw.githubusercontent.com/okadalabipr/biomass/master/resources/images/sensitivity_PcFos.png)

Control coefficients for integrated pc-Fos are shown by bars (blue, EGF; red, HRG). Numbers above bars indicate the reaction indices, and error bars correspond to simulation standard deviation.

## Citation

When using BioMASS, please cite:

- Imoto, H., Zhang, S. & Okada, M. A Computational Framework for Prediction and Analysis of Cancer Signaling Dynamics from RNA Sequencing Data—Application to the ErbB Receptor Signaling Pathway. _Cancers_. **12**, 2878 (2020). https://doi.org/10.3390/cancers12102878

  ```
  @article{imoto2020computational,
    title={A Computational Framework for Prediction and Analysis of Cancer Signaling Dynamics from RNA Sequencing Data—Application to the ErbB Receptor Signaling Pathway},
    author={Imoto, Hiroaki and Zhang, Suxiang and Okada, Mariko},
    journal={Cancers},
    volume={12},
    number={10},
    pages={2878},
    year={2020},
    publisher={Multidisciplinary Digital Publishing Institute}
  }
  ```

## Author

[Hiroaki Imoto](https://github.com/himoto)
