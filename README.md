# BioMASS

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Actions Status](https://github.com/okadalabipr/biomass/workflows/Tests/badge.svg)](https://github.com/okadalabipr/biomass/actions)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/okadalabipr/biomass.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/okadalabipr/biomass/context:python)
[![PyPI version](https://img.shields.io/pypi/v/biomass.svg?logo=PyPI&color=blue)](https://pypi.python.org/pypi/biomass/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/biomass.svg)](https://pypi.python.org/pypi/biomass/)

## Modeling and Analysis of Signaling Systems

<img align="left" src="https://github.com/okadalabipr/biomass_docs.jl/blob/master/docs/src/assets/logo.png?raw=true" width="300">

Mathematical modeling is a powerful method for the analysis of complex biological systems. Although there are many researches devoted on producing models to describe dynamical cellular signaling systems, most of these models are limited and do not cover multiple pathways. Therefore, there is a challenge to combine these models to enable understanding at a larger scale. Nevertheless, larger network means that it gets more difficult to estimate parameters to reproduce dynamic experimental data needed for deeper understanding of a system.

To overcome this problem, we developed BioMASS, a modeling platform tailored to optimizing mathematical models of biological processes. By using BioMASS, users can efficiently optimize kinetic parameters to fit user-defined models to experimental data, while performing analysis on reaction networks to predict critical components affecting cellular output.

## Features

BioMASS supports:

- Parameter Estimation of ODE Models
- Sensitivity Analysis
- Effective Visualization of Simulation Results

currently implimented for modeling immediate-early gene response ([Nakakuki _et al._, **_Cell_**, 2010](https://doi.org/10.1016/j.cell.2010.03.054)).

## Installation

The BioMASS library is available on [PyPI](https://pypi.org/project/biomass/).

```
$ pip install biomass
```

BioMASS supports Python 3.7 or newer.

## Create an executable model

```python
from biomass.models import Nakakuki_Cell_2010

Nakakuki_Cell_2010.show_properties()
```

```
Model properties
----------------
36 species
115 parameters, of which 75 to be estimated
```

```
model = Nakakuki_cell_2010.create()
```

## Parameter Estimation of ODE Models (_n_ = 1, 2, 3, · · ·)

Parameters are adjusted to minimize the distance between model simulation and experimental data.

```python
from biomass import optimize

optimize(
    model=model, start=1, options={
        "popsize": 3,
        "max_generation": 1000,
        "allowable_error": 0.5,
        "local_search_method": "DE",
    }
)
```

The temporary result will be saved in `out/n/` after each iteration.

Progress list: `out/n/optimization.log`

```
Generation1: Best Fitness = 1.726069e+00
Generation2: Best Fitness = 1.726069e+00
Generation3: Best Fitness = 1.726069e+00
Generation4: Best Fitness = 1.645414e+00
Generation5: Best Fitness = 1.645414e+00
Generation6: Best Fitness = 1.645414e+00
Generation7: Best Fitness = 1.645414e+00
Generation8: Best Fitness = 1.645414e+00
Generation9: Best Fitness = 1.645414e+00
Generation10: Best Fitness = 1.645414e+00
Generation11: Best Fitness = 1.645414e+00
Generation12: Best Fitness = 1.645414e+00
Generation13: Best Fitness = 1.645414e+00
Generation14: Best Fitness = 1.645414e+00
Generation15: Best Fitness = 1.645414e+00
Generation16: Best Fitness = 1.249036e+00
Generation17: Best Fitness = 1.171606e+00
Generation18: Best Fitness = 1.171606e+00
Generation19: Best Fitness = 1.171606e+00
Generation20: Best Fitness = 1.171606e+00
```

- If you want to continue from where you stopped in the last parameter search,

```python
from biomass import optimize_continue

optimize_continue(
    model=model, start=1, options={
        "popsize": 3,
        "max_generation": 1000,
        "allowable_error": 0.5,
        "local_search_method": "DE",
    }
)
```

- If you want to search multiple parameter sets (e.g., from 1 to 10) simultaneously,

```python
from biomass import optimize

optimize(
    model=model, start=1, end=10, options={
        "popsize": 5,
        "max_generation": 2000,
        "allowable_error": 0.5,
        "local_search_method": "mutation",
        "n_children": 50
    }
)
```

- Exporting optimized parameters in CSV format

```python
from biomass.result import OptimizationResults

res = OptimizationResults(model)
res.to_csv()
```

## Visualization of Simulation Results

```python
from biomass import run_simulation

run_simulation(model, viz_type='average', show_all=False, stdev=True)
```

**viz_type** : str

- `'average'`
  : The average of simulation results with parameter sets in `out/`.

- `'best'`
  : The best simulation result in `out/`, simulation with `best_fit_param`.

- `'original'`
  : Simulation with the default parameters and initial values defined in `set_model.py`.

- `'n(=1,2,...)'`
  : Use the parameter set in `out/n/`.
- `'experiment'`
  : Draw the experimental data written in `observable.py` without simulation results.

**show_all** : bool

- Whether to show all simulation results.

**stdev** : bool

- If True, the standard deviation of simulated values will be shown (only when `viz_type == 'average'`).

![simulation_average](https://github.com/okadalabipr/biomass_docs.jl/blob/master/docs/src/assets/simulation_average.png?raw=true)

Points (blue diamonds, EGF; red squares, HRG) denote experimental data, solid lines denote simulations

## Sensitivity Analysis

The single parameter sensitivity of each reaction is defined by<br>

_s<sub>i</sub>_(_q_(**v**),_v<sub>i</sub>_) = _∂_ ln(_q_(**v**)) / _∂_ ln(_v<sub>i</sub>_) = _∂_ _q_(**v**) / _∂_ _v<sub>i</sub>_ · _v<sub>i</sub>_ / _q_(**v**)

where _v<sub>i</sub>_ is the _i_<sup>th</sup> reaction rate, **v** is reaction vector **v** = (_v<sub>1</sub>_, _v<sub>2</sub>_, ...) and _q_(**v**) is a target function, e.g., time-integrated response, duration. Sensitivity coefficients were calculated using finite difference approximations with 1% changes in the reaction rates.

```python
from biomass import run_analysis

run_analysis(model, target='reaction', metric='integral', style='barplot')
```

**target** : str

- `'reaction'`
- `'initial_condition'`
- `'parameter'`

**metric** : str

- `'maximum'`
  : The maximum value.
- `'minimum'`
  : The minimum value.
- `'duration'`
  : The time it takes to decline below 10% of its maximum.
- `'integral'`
  : The integral of concentration over the observation time.

**style** : str

- `'barplot'`
- `'heatmap'`

![sensitivity_PcFos](https://github.com/okadalabipr/biomass_docs.jl/blob/master/docs/src/assets/sensitivity_PcFos.png?raw=true)

Control coefficients for integrated pc-Fos are shown by bars (blue, EGF; red, HRG). Numbers above bars indicate the reaction indices, and error bars correspond to simulation standard deviation.

## Citation

When using BioMASS, please cite:

- Imoto, H., Zhang, S. & Okada, M. A Computational Framework for Prediction and Analysis of Cancer Signaling Dynamics from RNA Sequencing Data—Application to the ErbB Receptor Signaling Pathway. _Cancers (Basel)_. **12**, 2878 (2020). https://doi.org/10.3390/cancers12102878
