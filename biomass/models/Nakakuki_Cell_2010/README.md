# Nakakuki_Cell_2010

This repository contains data and modeling code for the following paper:

Nakakuki, T. _et al._ Ligand-specific c-Fos expression emerges from the spatiotemporal control of ErbB network dynamics. _Cell_ **141**, 884–896 (2010). https://doi.org/10.1016/j.cell.2010.03.054

## BioModels

- [Nakakuki2010_CellFateDecision_Mechanistic](https://www.ebi.ac.uk/biomodels/BIOMD0000000250)

## Usage

1. Model construction

   ```python
   from biomass import create_model
   from biomass.models import Nakakuki_Cell_2010

   model = create_model(Nakakuki_Cell_2010.__package__, show_info=True)
   ```

   ↓

   ```
   Nakakuki_Cell_2010 information
   ------------------------------
   36 species
   115 parameters, of which 75 to be estimated
   ```

1. Parameter estimation

   ```python
   from tqdm import tqdm
   from biomass import optimize

   # Estimate 10 parameter sets
   for x_id in tqdm(range(1, 11)):
      optimize(model, x_id=x_id, optimizer_options={"workers": -1})
   ```

1. Visualization of simulation results

   ```python
   from biomass import run_simulation

   run_simulation(model, viz_type="average", show_all=False, stdev=True)
   ```

   ![simulation_average](https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/simulation_average.png)

1. Sensitivity analysis

   ```python
   from biomass import run_analysis

   run_analysis(model, target="reaction", metric='integral', style='barplot')
   ```

   ![sensitivity_PcFos](https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/sensitivity_PcFos.png)
