import os
import sys
from setuptools import setup, find_packages


def read_requirements():
    """Parse requirements from requirements.txt."""
    reqs_path = os.path.join(".", "requirements.txt")
    with open(reqs_path, "r") as f:
        requirements = [line.rstrip() for line in f]
    return requirements


def main():
    # Python version check.
    if sys.version_info[:2] < (3, 7):
        sys.exit("BioMASS requires at least Python version 3.7")

    # read version from file
    __version__ = ""
    version_path = os.path.join("biomass", "version.py")
    exec(open(version_path).read())

    # set long_description
    readme_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "README.md")
    with open(readme_path, "r", encoding="utf-8") as f:
        long_description = f.read()

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
        packages=find_packages(exclude=["docs", "tests"]),
        install_requires=read_requirements(),
        python_requires=">=3.7",
        classifiers=[
            "Intended Audience :: Science/Research",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    )


if __name__ == "__main__":
    main()