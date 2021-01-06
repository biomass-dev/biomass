# Lucarelli_CellSyst_2018

## Mathematical modeling of the TGF-β/Smad signaling pathway

Lucarelli, P. _et al._ Resolving the Combinatorial Complexity of Smad Protein Complex Formation and Its Link to Gene Expression. _Cell Syst._ **6**, 75-89.e11 (2018). https://doi.org/10.1016/j.cels.2017.11.010

## Simulate model

```python
from biomass.models import tgfb_smad
from biomass import run_simulation

model = tgfb_smad.create()

run_simulation(model, viz_type='original')
# Result : biomass/models/tgfb_smad/figure/simulation/original
```
