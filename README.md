# BioMASS

[![Actions Status](https://github.com/biomass-dev/biomass/workflows/Tests/badge.svg)](https://github.com/biomass-dev/biomass/actions)
[![Documentation Status](https://img.shields.io/readthedocs/biomass-core/latest.svg?logo=read%20the%20docs&logoColor=white&&label=Docs&version=latest)](https://biomass-core.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://img.shields.io/pypi/v/biomass.svg?logo=PyPI&logoColor=white)](https://pypi.python.org/pypi/biomass/)
[![License](https://img.shields.io/badge/License-Apache%202.0-green.svg?logo=apache)](https://opensource.org/licenses/Apache-2.0)
[![Downloads](https://pepy.tech/badge/biomass)](https://pepy.tech/project/biomass)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/biomass.svg?logo=Python&logoColor=white)](https://pypi.python.org/pypi/biomass/)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/biomass-dev/biomass.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/biomass-dev/biomass/context:python)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

<img align="left" src="https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/logo.png" width="300">

Mathematical modeling is a powerful method for the analysis of complex biological systems. Although there are many researches devoted on producing models to describe dynamical cellular signaling systems, most of these models are limited and do not cover multiple pathways. Therefore, there is a challenge to combine these models to enable understanding at a larger scale. Nevertheless, larger network means that it gets more difficult to estimate parameters to reproduce dynamic experimental data needed for deeper understanding of a system.

To overcome this problem, we developed _BioMASS_, a Python framework for **M**odeling and **A**nalysis of **S**ignaling **S**ystems. The BioMASS framework allows efficient optimization of multiple parameter sets simultaneously and generates the multiple parameter candidates that explain the signaling dynamics of interest. These parameter candidates can be further evaluated by their distribution and sensitivity analysis as a part of alternative information about the hidden regulatory mechanism of the system.

## Features

- Parameter estimation of ODE models
- Local sensitivity analysis
- Effective visualization of simulation results

## Documentation

Online documentation is available at https://biomass-core.readthedocs.io/.

## Installation

The BioMASS library is available at the [Python Package Index (PyPI)](https://pypi.org/project/biomass/).

```bash
$ pip install biomass
```

BioMASS supports Python 3.7 or newer.

## Example

### Parameter estimation

```python
from biomass import Model, optimize
from biomass.models import Nakakuki_Cell_2010

model = Model(Nakakuki_Cell_2010.__package__).create()

optimize(model, start=1, end=10)
```

![estimated_parameter_sets](https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/estimated_parameter_sets.png)

```python
from biomass import run_simulation

run_simulation(model, viz_type="average", stdev=True)
```

![simulation_average](https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/simulation_average.png)
Points (blue diamonds, EGF; red squares, HRG) denote experimental data, solid lines denote simulations.

### Sensitivity analysis

```python
from biomass import run_analysis

run_analysis(model, target="reaction", metric="integral", style="barplot")
```

![sensitivity_PcFos](https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/sensitivity_PcFos.png)

Control coefficients for integrated pc-Fos are shown by bars (blue, EGF; red, HRG). Numbers above bars indicate the reaction indices, and error bars correspond to simulation standard deviation.

## Citation

When using BioMASS, please cite the following paper:

- Imoto, H., Zhang, S. & Okada, M. A Computational Framework for Prediction and Analysis of Cancer Signaling Dynamics from RNA Sequencing Data—Application to the ErbB Receptor Signaling Pathway. _Cancers_ **12**, 2878 (2020). https://doi.org/10.3390/cancers12102878


## Author

[Hiroaki Imoto](https://github.com/himoto)

## License

[Apache License 2.0](https://github.com/biomass-dev/biomass/blob/master/LICENSE)
