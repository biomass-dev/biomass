import os
import shutil

from biomass import Text2Model


def test_preprocessing():
    if os.path.isdir(
        os.path.join(
            os.path.dirname(__file__),
            "text_files",
            "abc",
        )
    ):
        shutil.rmtree(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                "abc",
            )
        )


def test_text2model():
    abc = Text2Model(
        os.path.join(
            os.path.dirname(__file__),
            "text_files",
            "abc.txt",
        ),
    )
    abc.convert()
    # test thermodynamic restrictions
    desired = [["1", "2", "3", "4"]]
    actual = abc.restrictions
    assert len(actual) == 1
    assert set(desired[0]) == set(actual[0])


def test_cleanup():
    assert os.path.isdir(
        os.path.join(
            os.path.dirname(__file__),
            "text_files",
            "abc",
        )
    )
    shutil.rmtree(
        os.path.join(
            os.path.dirname(__file__),
            "text_files",
            "abc",
        )
    )
