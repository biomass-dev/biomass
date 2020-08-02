# Kholodenko_EurJBiochem_2000

## The quantitative computational model of the MAPK cascade
Kholodenko, B. N. Negative feedback and ultrasensitivity can bring about oscillations in the mitogen-activted protein kinase cascades. *Eur. J. Biochem* **267**, 1583–1588 (2000). https://doi.org/10.1046/j.1432-1327.2000.01197.x

## BioModels
- [Kholodenko2000 - Ultrasensitivity and negative feedback bring oscillations in MAPK cascade](https://www.ebi.ac.uk/biomodels/BIOMD0000000010)

## Simulate model

```python
from biomass.models import mapk_cascade
from biomass import run_simulation

run_simulation(mapk_cascade, viz_type='original')
# Result : biomass/models/mapk_cascade/figure/simulation/original
```