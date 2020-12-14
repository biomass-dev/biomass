import os
import sys
from setuptools import setup, find_packages


def read_file(file_path: str) -> str:
    """Read a file."""
    return open(file_path).read()


def main():
    # Python version check.
    if sys.version_info[:2] < (3, 7):
        sys.exit("BioMASS requires at least Python version 3.7")

    # read version from file
    __version__ = ""
    version_info = os.path.join("biomass", "version.py")
    exec(read_file(version_info))

    # set long_description and requirements
    here = os.path.abspath(os.path.dirname(__file__))
    long_description = read_file(os.path.join(here, "README.md"))
    requirements = read_file(os.path.join(here, "requirements.txt"))

    setup(
        name="biomass",
        version=__version__,
        description="A Python Framework for Modeling and Analysis of Signaling Systems",
        long_description=long_description,
        long_description_content_type="text/markdown",
        license="MIT",
        author="Hiroaki Imoto",
        author_email="himoto@protein.osaka-u.ac.jp",
        url="https://github.com/okadalabipr/biomass",
        packages=find_packages(exclude=["tests"]),
        install_requires=requirements.split(),
        python_requires=">=3.7",
        classifiers=[
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    )


if __name__ == "__main__":
    main()