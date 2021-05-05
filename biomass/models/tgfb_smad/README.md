# Lucarelli_CellSyst_2018

## Mathematical modeling of the TGF-β/Smad signaling pathway

Lucarelli, P. _et al._ Resolving the Combinatorial Complexity of Smad Protein Complex Formation and Its Link to Gene Expression. _Cell Syst._ **6**, 75-89.e11 (2018). https://doi.org/10.1016/j.cels.2017.11.010

## Simulation using BioMASS

Run this with Jupyter Notebook

```python
import os
from IPython.display import Image, display_png

from biomass.models import tgfb_smad
from biomass import Model, run_simulation

model = Model(tgfb_smad.__package__).create()

run_simulation(model, viz_type="original", save_format="png")

for observable in model.obs:
    with open(
        os.path.join(
            model.path,
            "figure",
            "simulation",
            "original",
            f"{observable}.png",
        ),
        mode="rb",
    ) as f:
        display_png(Image(f.read()))
```
