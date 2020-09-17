# BioMASS documentation

[![Actions Status](https://github.com/okadalabipr/biomass/workflows/Tests/badge.svg)](https://github.com/okadalabipr/biomass/actions)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/okadalabipr/biomass.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/okadalabipr/biomass/context:python)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

BioMASS is a user-friendly simulation tool for experimental biologists and currently implements the model of Nakakuki et al. 
BioMASS supports

- parameter estimation of ODE models
- sensitivity analysis
- effective visualization of simulation results

After model fitting, users can perform sensitivity analysis to identify critical parameters, species or regulations in the system of interest. 
By assessing parameter fitting and sensitivity analysis together, it is possible to identify critical or robust regions in the network, or the origin of heterogeneity generated from the network.

```@contents
    Pages = [
          "introduction/set_up.md",
          "introduction/overview.md",
          "introduction/import_model.md",
          "introduction/parameter_estimation.md",
          "introduction/visualization.md",
          "introduction/sensitivity_analysis.md",
    ]
    Depth = 3
```