import os
import shutil

from biomass import Text2Model

file_dir = os.path.join(os.path.dirname(__file__), "text_files")
txt_files = [name for name in os.listdir(file_dir) if "__" not in name]

skipped_files = ["duplicate_binding.txt", "typo_1.txt", "typo_2.txt", "unregistered_rule.txt"]


def test_preprocessing():
    for model in txt_files:
        model = model.split(".")[0]
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


def test_graph():
    for model_file in txt_files:
        if model_file in skipped_files:
            continue
        model = Text2Model(os.path.join(os.path.dirname(__file__), "text_files", model_file))
        model_path = os.path.join(file_dir, model_file.split(".")[0])
        model.register_word({"dissociate": ["releases"]})
        model.convert(overwrite=True)
        model.static_plot(save_dir=model_path, file_name="test.png")
        model.dynamic_plot(save_dir=model_path, file_name="test.html", show=False)
        assert os.stat(os.path.join(model_path, "test.png")).st_size > 1024 * 10
        assert os.stat(os.path.join(model_path, "test.html")).st_size > 1024 * 2


def test_cleanup():
    for model_file in txt_files:
        if model_file in skipped_files:
            continue
        model = model_file.split(".")[0]
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
