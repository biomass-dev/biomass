# Nakakuki_Cell_2010
This repository contains data and modeling code for the following paper:

Nakakuki, T. *et al.* Ligand-specific c-Fos expression emerges from the spatiotemporal control of ErbB network dynamics. *Cell* **141**, 884–896 (2010). https://doi.org/10.1016/j.cell.2010.03.054

## BioModels
- [Nakakuki2010_CellFateDecision_Mechanistic](https://www.ebi.ac.uk/biomodels/BIOMD0000000250)

## Usage
1. Parameter estimation
    ```python
    from biomass.models import Nakakuki_Cell_2010
    from biomass import optimize

    # Estimate 10 parameter sets simultaneously
    optimize(Nakakuki_Cell_2010, 1, 10, max_generation=10000, allowable_error=0.5)
    ```

1. Visualization of simulation results
    ```python
    from biomass.models import Nakakuki_Cell_2010
    from biomass import run_simulation

    run_simulation(Nakakuki_Cell_2010, viz_type='average', show_all=False, stdev=True)
    ```

1. Sensitivity analysis
    ```python
    from biomass.models import Nakakuki_Cell_2010
    from biomass import run_analysis
    
    run_analysis(Nakakuki_Cell_2010, target='reaction', metric='integral', style='barplot')
    ```