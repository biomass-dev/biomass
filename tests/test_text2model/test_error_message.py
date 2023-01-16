import os
import shutil

import pytest

from biomass import Text2Model
from biomass.construction.reaction_rules import DetectionError


def test_preprocessing():
    for i in ["1", "2"]:
        if os.path.isdir(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                f"typo_{i}",
            )
        ):
            shutil.rmtree(
                os.path.join(
                    os.path.dirname(__file__),
                    "text_files",
                    f"typo_{i}",
                )
            )


def test_typo_detection():
    for i in ["1", "2"]:
        typo = Text2Model(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                f"typo_{i}.txt",
            ),
        )
        with pytest.raises(DetectionError) as e:
            typo.convert()
        assert "Maybe: 'binds'." in str(e.value)


def test_unregistered_rule():
    ubiquitination = Text2Model(
        os.path.join(
            os.path.dirname(__file__),
            "text_files",
            "unregistered_rule.txt",
        )
    )
    with pytest.raises(DetectionError) as e:
        ubiquitination.convert()
    assert "Unregistered words in line1: A ubiquitinates B --> uB" in str(e.value)


def test_cleanup():
    for i in ["1", "2"]:
        assert os.path.isdir(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                f"typo_{i}",
            )
        )
        shutil.rmtree(
            os.path.join(
                os.path.dirname(__file__),
                "text_files",
                f"typo_{i}",
            )
        )
    assert os.path.isdir(
        os.path.join(
            os.path.dirname(__file__),
            "text_files",
            "unregistered_rule",
        )
    )
    shutil.rmtree(
        os.path.join(
            os.path.dirname(__file__),
            "text_files",
            "unregistered_rule",
        )
    )
