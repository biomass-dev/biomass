import os
import sys
from setuptools import setup, find_packages


def get_version() -> str:
    """Read version from file"""
    version_filepath = os.path.join(os.path.dirname(__file__), "biomass", "version.py")
    with open(version_filepath) as f:
        for line in f:
            if line.startswith("__version__"):
                return line.strip().split()[-1][1:-1]


def main():
    # Python version check.
    if sys.version_info[:2] < (3, 7):
        sys.exit("BioMASS requires at least Python version 3.7")

    # set long_description and requirements
    here = os.path.abspath(os.path.dirname(__file__))
    long_description = open(os.path.join(here, "README.md")).read()
    requirements = open(os.path.join(here, "requirements.txt")).read()

    setup(
        name="biomass",
        version=get_version(),
        description="A Python Framework for Modeling and Analysis of Signaling Systems",
        long_description=long_description,
        long_description_content_type="text/markdown",
        license="MIT",
        author="Hiroaki Imoto",
        author_email="himoto@protein.osaka-u.ac.jp",
        url="https://github.com/okadalabipr/biomass",
        packages=find_packages(exclude=["tests"]),
        install_requires=requirements.splitlines(),
        python_requires=">=3.7",
        keywords=["systems", "biology", "modeling", "optimization", "sensitivity", "analysis"],
        classifiers=[
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
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