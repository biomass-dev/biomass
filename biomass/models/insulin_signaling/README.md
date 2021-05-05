# Kubota_MolCell_2012

## The Insulin-Dependent AKT Pathway Model

Kubota, H. _et al._ Temporal Coding of Insulin Action through Multiplexing of the AKT Pathway. _Mol. Cell_ **46**, 820–832 (2012). https://doi.org/10.1016/j.molcel.2012.04.018

## BioModels

- [Kubota2012_InsulinAction_AKTpathway](https://www.ebi.ac.uk/biomodels/MODEL1204060000)

## Simulation using BioMASS

Run this with Jupyter Notebook

```python
import os
from IPython.display import Image, display_png

from biomass.models import insulin_signaling
from biomass import Model, run_simulation

model = Model(insulin_signaling.__package__).create()

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
