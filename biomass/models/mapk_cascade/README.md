# Kholodenko_EurJBiochem_2000

## The quantitative computational model of the MAPK cascade

Kholodenko, B. N. Negative feedback and ultrasensitivity can bring about oscillations in the mitogen-activted protein kinase cascades. _Eur. J. Biochem_ **267**, 1583–1588 (2000). https://doi.org/10.1046/j.1432-1327.2000.01197.x

## BioModels

- [Kholodenko2000 - Ultrasensitivity and negative feedback bring oscillations in MAPK cascade](https://www.ebi.ac.uk/biomodels/BIOMD0000000010)

## Simulation using BioMASSS

Run this with Jupyter Notebook

```python
import os
from IPython.display import Image, display_png

from biomass.models import mapk_cascade
from biomass import Model, run_simulation

model = Model(mapk_cascade.__package__).create()

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
