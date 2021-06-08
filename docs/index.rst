|Actions Status| |Documentation Status| |PyPI version| |License| |Downloads| |Python versions| |Code quality| |Code style|

=====================
BioMASS documentation
=====================

Mathematical modeling is a powerful method for the analysis of complex biological systems.
Although there are many researches devoted on producing models to describe dynamical cellular signaling systems, most of these models are limited and do not cover multiple pathways.
Therefore, there is a challenge to combine these models to enable understanding at a larger scale.
Nevertheless, larger network means that it gets more difficult to estimate parameters to reproduce dynamic experimental data needed for deeper understanding of a system.

To overcome this problem, we developed *BioMASS*, a Python framework for Modeling and Analysis of Signaling Systems.
The BioMASS framework allows efficient optimization of multiple parameter sets simultaneously and generates the multiple parameter candidates that explain the signaling dynamics of interest.
These parameter candidates can be further evaluated by their distribution and sensitivity analysis as a part of alternative information about the hidden regulatory mechanism of the system.

Contents
========

.. toctree::
   :maxdepth: 2

   about
   installation
   models
   usage
   modules/index
   citation

.. |Actions Status| image:: https://github.com/biomass-dev/biomass/workflows/Tests/badge.svg
   :target: https://github.com/biomass-dev/biomass/actions
   :alt: Actions Status

.. |Documentation Status| image:: https://img.shields.io/readthedocs/biomass-core/latest.svg?logo=read%20the%20docs&logoColor=white&&label=Docs&version=latest
   :target: https://biomass-core.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |PyPI version| image:: https://img.shields.io/pypi/v/biomass.svg?logo=PyPI&logoColor=white
   :target: https://pypi.python.org/pypi/biomass/
   :alt: PyPI version

.. |License| image:: https://img.shields.io/badge/License-Apache%202.0-green.svg
   :target: https://opensource.org/licenses/Apache-2.0
   :alt: License

.. |Downloads| image:: https://pepy.tech/badge/biomass
   :target: https://pepy.tech/project/biomass
   :alt: Downloads

.. |Python versions| image:: https://img.shields.io/pypi/pyversions/biomass.svg?logo=Python&logoColor=white
   :target: https://pypi.python.org/pypi/biomass/
   :alt: Python versions

.. |Code quality| image:: https://img.shields.io/lgtm/grade/python/g/biomass-dev/biomass.svg?logo=lgtm&logoWidth=18
   :target: https://lgtm.com/projects/g/biomass-dev/biomass/context:python
   :alt: Code quality: Python

.. |Code style| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Code style: black