import os
import shutil
from typing import Callable

import numpy as np
import pytest

from biomass import Model, Text2Model, run_simulation
from biomass.construction.thermodynamic_restrictions import DuplicateError


def test_preprocessing():
    for model in [
        "michaelis_menten",
        "Kholodenko1999",
        "michaelis_menten_jl",
        "Kholodenko1999_jl",
    ]:
        if os.path.isdir(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                model,
            )
        ):
            shutil.rmtree(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    model,
                )
            )


def test_text2model():
    for model in ["michaelis_menten", "Kholodenko1999"]:
        for lang in ["python", "julia"]:
            if model == "michaelis_menten":
                mm_kinetics = Text2Model(
                    os.path.join(
                        os.path.dirname(__file__),
                        "text_files",
                        f"{model}.txt",
                    ),
                    similarity_threshold=0.7,
                    lang=lang,
                )
                mm_kinetics.register_word({"dissociate": ["releases"]})
                mm_kinetics.convert()
            elif model == "Kholodenko1999":
                mapk_cascade = Text2Model(
                    os.path.join(
                        os.path.dirname(__file__),
                        "text_files",
                        f"{model}.txt",
                    ),
                    lang=lang,
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

    with pytest.raises(DuplicateError):
        Text2Model(
            os.path.join(os.path.dirname(__file__), "text_files", "duplicate_binding.txt")
        ).convert()


def test_run_simulation():
    try:
        from .text_files import Kholodenko1999, michaelis_menten

        _packaging: Callable[[str], str] = lambda model_name: ".".join(
            ["tests.test_text2model.text_files", model_name]
        )
        for model_name in ["Kholodenko1999", "michaelis_menten"]:
            model = Model(_packaging(model_name)).create()
            run_simulation(model)
            simulated_values = np.load(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    model_name,
                    "simulation_data",
                    "simulations_original.npy",
                )
            )
            assert np.isfinite(simulated_values).all()
            if model_name == "Kholodenko1999":
                expected = np.load(
                    os.path.join(
                        os.path.dirname(__file__),
                        "simulations_original_BN.npy",
                    )
                )
                assert np.allclose(simulated_values, expected)
    except ImportError as e:
        print(e)


def test_text2markdown():
    for model in ["michaelis_menten", "Kholodenko1999"]:
        if model == "michaelis_menten":
            mm_kinetics = Text2Model(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    f"{model}.txt",
                )
            )
            mm_kinetics.register_word({"dissociate": ["releases"]})
            mm_kinetics.to_markdown(num_reactions=2)
        elif model == "Kholodenko1999":
            mapk_cascade = Text2Model(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    f"{model}.txt",
                )
            )
            mapk_cascade.to_markdown(num_reactions=25)
        assert os.path.isfile(os.path.join("markdown", model, "rate.md"))
        assert os.path.isfile(os.path.join("markdown", model, "diffeq.md"))
        shutil.rmtree("markdown")


def test_julia_models():
    necessities = [
        os.path.join("name2idx", "parameters.jl"),
        os.path.join("name2idx", "species.jl"),
        "ode.jl",
        "observable.jl",
        "simulation.jl",
        "experimental_data.jl",
        "search_param.jl",
        "problem.jl",
    ]
    for model in ["michaelis_menten", "Kholodenko1999"]:
        for file in necessities:
            assert os.path.isfile(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    f"{model}_jl",
                    file,
                )
            )


def test_cleanup():
    for model in [
        "michaelis_menten",
        "Kholodenko1999",
        "michaelis_menten_jl",
        "Kholodenko1999_jl",
        "duplicate_binding",
    ]:
        assert os.path.isdir(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                model,
            )
        )
        shutil.rmtree(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                model,
            )
        )
