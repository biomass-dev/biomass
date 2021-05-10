import os
import shutil
from distutils.dir_util import copy_tree

import pytest

from biomass import Model, run_simulation
from biomass.models import Nakakuki_Cell_2010

# from biomass import run_analysis


os.makedirs(
    os.path.join(
        "biomass",
        "models",
        "Nakakuki_Cell_2010",
        "out",
    ),
    exist_ok=True,
)
copy_tree(
    os.path.join(
        "tests",
        "out",
    ),
    os.path.join(
        "biomass",
        "models",
        "Nakakuki_Cell_2010",
        "out",
    ),
)

model = Model(Nakakuki_Cell_2010.__package__).create()

for dir in ["figure", "simulation_data", "sensitivity_coefficients"]:
    if os.path.isdir(os.path.join(model.path, dir)):
        shutil.rmtree(os.path.join(model.path, dir))


def test_run_simulation():
    with pytest.warns(RuntimeWarning) as record:
        run_simulation(model, viz_type="original")
    # check that the message matches
    assert record[-1].message.args[0] in [
        "Simulation failed. #original",
        "invalid value encountered in true_divide",
    ]
    run_simulation(model, viz_type="average", stdev=True)
    assert os.path.isfile(
        os.path.join(
            model.path,
            "simulation_data",
            "simulations_all.npy",
        )
    )


def test_cleanup():
    for dir in ["figure", "simulation_data", "out", "sensitivity_coefficients"]:
        if os.path.isdir(os.path.join(model.path, dir)):
            shutil.rmtree(os.path.join(model.path, dir))
