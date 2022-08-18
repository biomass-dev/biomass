import os
import shutil
from typing import Callable

import numpy as np

from biomass import Text2Model, create_model, run_simulation


def test_preprocessing():
    if os.path.isdir(
        os.path.join(
            os.path.dirname(__file__),
            "text_files",
            "Kholodenko1999_2",
        )
    ):
        shutil.rmtree(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                "Kholodenko1999_2",
            )
        )


def test_text2model():
    mapk_cascade = Text2Model(
        os.path.join(
            os.path.dirname(__file__),
            "text_files",
            "Kholodenko1999_2.txt",
        ),
    )
    mapk_cascade.convert()
    # test thermodynamic restrictions
    desired = [
        {"10", "11", "12", "9"},
        {"15", "17", "18", "21"},
        {"18", "19", "20", "22"},
        {"18", "19", "20", "22"},
        {"12", "17", "19", "24"},
        {"15", "20", "23", "24"},
        {"12", "21", "22", "23"},
    ]
    actual = [set(restriction) for restriction in mapk_cascade.restrictions]
    for l1 in actual:
        assert l1 in desired
    for l2 in desired:
        assert l2 in actual


def test_run_simulation():
    try:
        from .text_files import Kholodenko1999_2

        _packaging: Callable[[str], str] = lambda model_name: ".".join(
            ["tests.test_text2model.text_files", model_name]
        )
        model = create_model(_packaging("Kholodenko1999_2"))
        run_simulation(model)
        simulated_values = np.load(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                "Kholodenko1999_2",
                "simulation_data",
                "simulations_original.npy",
            )
        )
        assert np.isfinite(simulated_values).all()
        expected = np.load(
            os.path.join(
                os.path.dirname(__file__),
                "simulations_original_BN.npy",
            )
        )
        assert np.allclose(simulated_values, expected)
    except ImportError as e:
        print(e)


def test_cleanup():
    if os.path.isdir(
        os.path.join(
            os.path.dirname(__file__),
            "text_files",
            "Kholodenko1999_2",
        )
    ):
        shutil.rmtree(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                "Kholodenko1999_2",
            )
        )
