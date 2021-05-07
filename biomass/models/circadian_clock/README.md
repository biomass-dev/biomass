# Leloup_PNAS_2003

## Computational Model for the Mammalian Circadian Clock

Leloup, J.-C. & Goldbeter, A. Toward a detailed computational model for the mammalian circadian clock. _Proc. Natl. Acad. Sci._ **100**, 7051–6 (2003). https://doi.org/10.1073/pnas.1132112100

## BioModels

- [Leloup2003_CircClock_DD](https://www.ebi.ac.uk/biomodels/BIOMD0000000073)

## Simulation using BioMASS

Run this with Jupyter Notebook

```python
import os
from IPython.display import Image, display_png

from biomass.models import circadian_clock
from biomass import Model, run_simulation

model = Model(circadian_clock.__package__).create()

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
