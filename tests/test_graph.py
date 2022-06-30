import os

import pytest

from biomass import create_model


def test_graph():
    model_lib_list = [name for name in os.listdir("biomass/models") if "__" not in name]
    for model_name in model_lib_list:
        if model_name in ["prolif_quies", "g1s_transition"]:
            model = create_model(".".join(("biomass", "models", model_name)))
            with pytest.warns(UserWarning):
                model.static_plot(os.path.join('biomass', 'models', model_name),os.path.join(model_name + "test.png"))
                assert (
                    os.stat(os.path.join(model.path, model_name + "test.png")).st_size > 1024 * 10
                )
                model.dynamic_plot(os.path.join('biomass', 'models', model_name),os.path.join(model_name + "test.html"), show=False)
                assert (
                    os.stat(os.path.join(model.path, model_name + "test.html")).st_size > 1024 * 2
                )
        else:
            model = create_model(".".join(("biomass", "models", model_name)))
            model.static_plot(os.path.join('biomass', 'models', model_name),os.path.join(model_name + "test.png"))
            model.dynamic_plot(os.path.join('biomass', 'models', model_name),os.path.join(model_name + "test.html"), show=False)
            assert os.stat(os.path.join(model.path, model_name + "test.png")).st_size > 1024 * 10
            assert os.stat(os.path.join(model.path, model_name + "test.html")).st_size > 1024 * 2
        os.remove(os.path.join(model.path, model_name + "test.png"))
        os.remove(os.path.join(model.path, model_name + "test.html"))