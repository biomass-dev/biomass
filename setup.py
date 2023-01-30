import os
import sys
from typing import List

from setuptools import find_packages, setup

try:
    from biomass import __author__, __email__, __maintainer__
except ImportError:
    __author__ = __maintainer__ = "Hiroaki Imoto"
    __email__ = "himoto@protein.osaka-u.ac.jp"


def get_version() -> str:

    version_filepath = os.path.join(os.path.dirname(__file__), "biomass", "version.py")
    with open(version_filepath) as f:
        for line in f:
            if line.startswith("__version__"):
                return line.strip().split()[-1][1:-1]
    assert False


def get_long_description() -> str:

    readme_filepath = os.path.join(os.path.dirname(__file__), "README.md")
    with open(readme_filepath) as f:
        return f.read()


def get_install_requires() -> List[str]:

    requirements_filepath = os.path.join(os.path.dirname(__file__), "requirements.txt")
    with open(requirements_filepath) as f:
        return f.read().splitlines()


def setup_package():
    # Python version check.
    if sys.version_info[:2] < (3, 8):
        raise RuntimeError("biomass requires at least Python version 3.8")

    setup(
        name="biomass",
        version=get_version(),
        description="A Python Framework for Modeling and Analysis of Signaling Systems",
        long_description=get_long_description(),
        long_description_content_type="text/markdown",
        license="Apache 2.0",
        author=__author__,
        author_email=__email__,
        maintainer=__maintainer__,
        maintainer_email=__email__,
        url="https://github.com/biomass-dev/biomass",
        download_url="https://github.com/biomass-dev/biomass/releases",
        project_urls={
            "Documentation": "https://biomass-core.readthedocs.io/en/latest/",
            "Source Code": "https://github.com/biomass-dev/biomass",
            "Bug Tracker": "https://github.com/biomass-dev/biomass/issues",
        },
        packages=find_packages(exclude=["tests", "docs"]),
        install_requires=get_install_requires(),
        extras_require={
            "graph": [
                "pygraphviz>=1.9",
                "pyvis>=0.2.1,<0.3",
            ],
            "dev": [
                "black==22.12.0",
                "flake8==6.0.0",
                "isort==5.12.0",
                "pre-commit",
                "pytest",
            ],
            "docs": [
                "sphinx>=1.7,<6",
                "sphinx_rtd_theme>=0.3",
                "sphinx_autodoc_typehints>=1.10",
                "sphinxcontrib-bibtex>=2.2",
            ],
        },
        python_requires=">=3.8",
        keywords=[
            "systems",
            "biology",
            "kinetic",
            "modeling",
            "biochemical",
            "reaction",
            "network",
            "simulation",
            "ode",
        ],
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: Apache Software License",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3 :: Only",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Software Development",
            "Topic :: Software Development :: Libraries",
            "Topic :: Software Development :: Libraries :: Python Modules",
        ],
    )


if __name__ == "__main__":
    setup_package()
