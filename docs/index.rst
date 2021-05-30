.. image:: https://github.com/okadalabipr/biomass/workflows/Tests/badge.svg
   :target: https://github.com/okadalabipr/biomass/actions

.. image:: https://img.shields.io/pypi/v/biomass.svg?logo=PyPI&logoColor=white
   :target: https://pypi.python.org/pypi/biomass/

.. image:: https://img.shields.io/badge/License-Apache%202.0-green.svg
   :target: https://opensource.org/licenses/Apache-2.0
   
.. image:: https://pepy.tech/badge/biomass
   :target: https://pepy.tech/project/biomass

.. image:: https://img.shields.io/pypi/pyversions/biomass.svg?logo=Python&logoColor=white
   :target: https://pypi.python.org/pypi/biomass/

.. image:: https://img.shields.io/lgtm/grade/python/g/okadalabipr/biomass.svg?logo=lgtm&logoWidth=18
   :target: https://lgtm.com/projects/g/okadalabipr/biomass/context:python

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black

BioMASS
=======

.. image:: https://raw.githubusercontent.com/okadalabipr/biomass/master/docs/_static/img/logo.png
   :width: 300px
   :align: left

Mathematical modeling is a powerful method for the analysis of complex biological systems.
Although there are many researches devoted on producing models to describe dynamical cellular signaling systems, most of these models are limited and do not cover multiple pathways.
Therefore, there is a challenge to combine these models to enable understanding at a larger scale.
Nevertheless, larger network means that it gets more difficult to estimate parameters to reproduce dynamic experimental data needed for deeper understanding of a system.

To overcome this problem, we developed *BioMASS*, a Python framework for Modeling and Analysis of Signaling Systems.
The BioMASS framework allows efficient optimization of multiple parameter sets simultaneously and generates the multiple parameter candidates that explain the signaling dynamics of interest.
These parameter candidates can be further evaluated by their distribution and sensitivity analysis as a part of alternative information about the hidden regulatory mechanism of the system.

.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   models
   usage
   modules/index
   citation

Author & Maintainer
-------------------

`Hiroaki Imoto`_

License
-------

`Apache License 2.0`_

.. _`Hiroaki Imoto`: https://github.com/himoto
.. _`Apache License 2.0`: https://github.com/okadalabipr/biomass/blob/master/LICENSE