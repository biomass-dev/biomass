import os
import sys
from typing import List

from setuptools import find_packages, setup


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


def main():
    # Python version check.
    if sys.version_info[:2] < (3, 7):
        sys.exit("biomass requires at least Python version 3.7")

    setup(
        name="biomass",
        version=get_version(),
        description="A Python Framework for Modeling and Analysis of Signaling Systems",
        long_description=get_long_description(),
        long_description_content_type="text/markdown",
        license="MIT",
        author="Hiroaki Imoto",
        author_email="himoto@protein.osaka-u.ac.jp",
        url="https://github.com/okadalabipr/biomass",
        download_url="https://pypi.org/project/biomass/",
        packages=find_packages(exclude=["tests"]),
        install_requires=get_install_requires(),
        extras_require={
            "dev": [
                "black==20.8b1",
                "flake8",
                "isort",
                "pre-commit",
                "pytest",
            ]
        },
        python_requires=">=3.7",
        keywords=[
            "systems",
            "biology",
            "modeling",
            "optimization",
            "sensitivity",
            "analysis",
        ],
        classifiers=[
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3 :: Only",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Software Development",
            "Topic :: Software Development :: Libraries",
            "Topic :: Software Development :: Libraries :: Python Modules",
        ],
    )


if __name__ == "__main__":
    main()
