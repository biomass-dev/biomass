import os
import shutil

import pytest

from biomass import create_model


def test_graph():
    model_lib_list = [name for name in os.listdir("biomass/models") if "__" not in name]
    for model_name in model_lib_list:
        if model_name in ["prolif_quies", "g1s_transition"]:
            model = create_model(".".join(("biomass", "models", model_name)))
            with pytest.warns(UserWarning):
                model.to_graph(os.path.join(model_name + "test.png"))
                assert (
                    os.stat(
                        os.path.join(model.path, model_name + "test.png")
                    ).st_size
                    > 1024 * 10
                )
            continue
        model = create_model(".".join(("biomass", "models", model_name)))
        model.to_graph(os.path.join(model_name + "test.png"))
        assert (
            os.stat(os.path.join(model.path, model_name + "test.png")).st_size
            > 1024 * 10
        )
        os.remove(os.path.join(model.path, model_name + "test.png"))
