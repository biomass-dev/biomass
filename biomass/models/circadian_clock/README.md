# Leloup_PNAS_2003

## Computational Model for the Mammalian Circadian Clock
Leloup, J.-C. & Goldbeter, A. Toward a detailed computational model for the mammalian circadian clock. *Proc. Natl. Acad. Sci.* **100**, 7051–6 (2003). https://doi.org/10.1073/pnas.1132112100

## BioModels
- [Leloup2003_CircClock_DD](https://www.ebi.ac.uk/biomodels/BIOMD0000000073)

## Simulate model

```python
from biomass.models import circadian_clock
from biomass import run_simulation

run_simulation(circadian_clock, viz_type='original')
# Result : biomass/models/circadian_clock/figure/simulation/original
```