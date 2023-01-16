# Kholodenko_EurJBiochem_2000

## The quantitative computational model of the MAPK cascade

Kholodenko, B. N. Negative feedback and ultrasensitivity can bring about oscillations in the mitogen-activted protein kinase cascades. _Eur. J. Biochem_ **267**, 1583–1588 (2000). https://doi.org/10.1046/j.1432-1327.2000.01197.x

## BioModels

- [Kholodenko2000 - Ultrasensitivity and negative feedback bring oscillations in MAPK cascade](https://www.ebi.ac.uk/biomodels/BIOMD0000000010)

## Run simulation using BioMASSS

```python
import os

from biomass import create_model, run_simulation
from biomass.models import copy_to_current

copy_to_current("mapk_cascade")
model = create_model("mapk_cascade")

run_simulation(model)
```

<img align="left" src="./mapk_cascade.png" width="400px">
