# Oppelt_npjSystBiolAppl_2018

## Dynamic pathway model of TNFα-induced NFκB signal transduction

Oppelt, A. _et al._ Model-based identification of TNFα-induced IKKβ-mediated and IκBα-mediated regulation of NFκB signal transduction as a tool to quantify the impact of drug-induced liver injury compounds. _npj Syst. Biol. Appl._ **4**, 23 (2018). https://doi.org/10.1038/s41540-018-0058-z

## Simulation using BioMASS

Run this with Jupyter Notebook

```python
import os
from IPython.display import Image, display_png

from biomass.models import nfkb_pathway
from biomass import Model, run_simulation

model = Model(nfkb_pathway.__package__).create()

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
