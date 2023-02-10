About
=====

*BioMASS* is a computational framework for modeling and analysis of biological signaling systems in Python.
It provides useful tools for model construction, numerical simulation, parameter estimation, network analysis, and result visualization.

Example
-------

Text file (``michaelis_menten.txt``):

.. code-block::
    :linenos:
    
    E + S ⇄ ES | kf=0.003, kr=0.001 | E=100, S=50
    ES → E + P | kf=0.002
    
    @obs Substrate: u[S]
    @obs E_free: u[E]
    @obs E_total: u[E] + u[ES]
    @obs Product: u[P]
    @obs Complex: u[ES]
    
    @sim tspan: [0, 100]

Text-to-model conversion:

.. code-block:: python

  >>> from biomass import Text2Model, create_model, run_simulation
  >>> description = Text2Model("michaelis_menten.txt")
  >>> description.convert()
  Model information
  -----------------
  2 reactions
  4 species
  4 parameters

  >>> model = create_model("michaelis_menten")
  >>> run_simulation(model)
    
Output:

.. image:: https://raw.githubusercontent.com/pasmopy/pasmopy/master/docs/_static/img/michaelis_menten_sim.png

For an advanced model, see `a mechanistic model of the c-Fos expression network dynamics <https://biomass-core.readthedocs.io/en/latest/tutorial/cfos.html>`_.

License
-------

The software is released under the `Apache License 2.0 <https://opensource.org/licenses/Apache-2.0>`_.
For details, see the `LICENSE <https://github.com/biomass-dev/biomass/blob/master/LICENSE>`_ file in the biomass repository.

Author
------

`Hiroaki Imoto <https://github.com/himoto>`_

Contact
-------

Please contact me with any questions or comments via `Issues`_ |  `Discussions`_ on GitHub.
You can also always send me an `email <mailto:himoto@protein.osaka-u.ac.jp>`_.

Any contributions to BioMASS are more than welcome!

.. _Issues: https://github.com/biomass-dev/biomass/issues
.. _Discussions: https://github.com/biomass-dev/biomass/discussions
