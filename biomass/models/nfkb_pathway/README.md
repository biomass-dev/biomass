# Oppelt_npjSystBiolAppl_2018

## Dynamic pathway model of TNFα-induced NFκB signal transduction

Oppelt, A. _et al._ Model-based identification of TNFα-induced IKKβ-mediated and IκBα-mediated regulation of NFκB signal transduction as a tool to quantify the impact of drug-induced liver injury compounds. _npj Syst. Biol. Appl._ **4**, 23 (2018). https://doi.org/10.1038/s41540-018-0058-z

## Simulate model

```python
from biomass.models import nfkb_pathway
from biomass import run_simulation

model = nfkb_pathway.create()

run_simulation(model, viz_type='original')
# Result : biomass/models/nfkb_pathway/figure/simulation/original
```
