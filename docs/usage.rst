Usage
=====

.. image:: https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/overview.png

Model preparation
-----------------

A brief description of each file/folder is below:

======================= ========================================================================================================
Name                    Content
======================= ========================================================================================================
``name2idx/``           Names of model parameters and species
``set_model.py``        Differential equation, parameters and initial condition
``observalbe.py``       Observables, simulations and experimental data
``viz.py``              Plotting parameters for customizing figure properties
``set_search_param.py`` Model parameters to optimize and search region
``fitness.py``          An objective function to be minimized, i.e., the distance between model simulation and experimental data
``reaction_network.py`` Reaction indices grouped according to biological processes
======================= ========================================================================================================

.. code-block:: shell

    $ git clone https://github.com/biomass-dev/biomass
    $ cd biomass
    $ python

.. code-block:: python

    >>> from biomass import create_model
    >>> from biomass.models import Nakakuki_Cell_2010
    >>> model = create_model(Nakakuki_Cell_2010.__package__, show_info=True)
    Nakakuki_Cell_2010 information
    ------------------------------
    36 species
    115 parameters, of which 75 to be estimated

.. note::
    `pasmopy.Text2Model <https://pasmopy.readthedocs.io/en/latest/model_development.html>`_ allows you to build a BioMASS model from text :cite:p:`imoto2022text`.
    You simply describe biochemical reactions and the molecular mechanisms extracted from text are converted into an executable model.

Parameter estimation
--------------------

Using :func:`~biomass.core.optimize` function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parameters are adjusted to minimize the distance between model simulation and experimental data.

* Set simulation conditions and the corresponding experimental data in ``observable.py``
* Define an objective function to be minimized (:func:`objective`) in ``fitness.py``
* Set lower/upper bounds of parameters to be estimated in ``set_search_param.py``

.. code-block:: python

    from biomass import optimize
    
    optimize(model, x_id=1, optimizer_options={"workers": -1})

The temporary result will be saved in ``out/_tmp{n}/`` after each iteration.

Progress list: ``out/_tmp{n}/optimization.log``::

    differential_evolution step 1: f(x)= 4.96181
    differential_evolution step 2: f(x)= 3.555
    differential_evolution step 3: f(x)= 2.50626
    differential_evolution step 4: f(x)= 2.00657
    differential_evolution step 5: f(x)= 1.83556
    differential_evolution step 6: f(x)= 1.28031
    differential_evolution step 7: f(x)= 0.973207
    differential_evolution step 8: f(x)= 0.741667
    differential_evolution step 9: f(x)= 0.741667
    differential_evolution step 10: f(x)= 0.735682
    differential_evolution step 11: f(x)= 0.717266
    differential_evolution step 12: f(x)= 0.603178
    differential_evolution step 13: f(x)= 0.56934
    differential_evolution step 14: f(x)= 0.56934
    differential_evolution step 15: f(x)= 0.549331
    differential_evolution step 16: f(x)= 0.459069
    differential_evolution step 17: f(x)= 0.447772
    differential_evolution step 18: f(x)= 0.430385
    differential_evolution step 19: f(x)= 0.37085
    differential_evolution step 20: f(x)= 0.37085

.. note::
    For detailed information about ``optimizer_options``, please refer to `scipy docs <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html>`_.

Data export and visualization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from biomass.result import OptimizationResults

    res = OptimizationResults(model)
    # Export estimated parameters in CSV format
    res.to_csv()
    # Visualize estimated parameter sets
    res.savefig(figsize=(16,5), boxplot_kws={"orient": "v"})

.. image:: https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/estimated_parameter_sets.png

.. code-block:: python

    # Visualize objective function traces for different optimization runs.
    res.trace_obj()

.. image:: https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/obj_func_trace.png

Visualization of simulation results
-----------------------------------

.. code-block:: python

    from biomass import run_simulation

    run_simulation(model, viz_type='average', show_all=False, stdev=True)

.. image:: https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/simulation_average.png

Points (blue diamonds, EGF; red squares, HRG) denote experimental data, solid lines denote simulations.

Sensitivity analysis
--------------------

Sensitivity analysis examines how perturbations to the processes in the model affect the quantity of interest, e.g., the integral of the pc-Fos concentration.

.. code-block:: python

    from biomass import run_analysis

    run_analysis(model, target='reaction', metric='integral', style='barplot', options={'overwrite': True})

The single parameter sensitivity of each reaction is defined by

.. math:: C^{M}_{i} = d \ln{M} / d \ln{v_{i}}

where v\ :sub:`i`\  is the i\ :sup:`th`\  reaction rate, v is reaction vector v = (v\ :sub:`1`\, v\ :sub:`2`\, ...) and M is a signaling metric, e.g., time-integrated response, duration.
Sensitivity coefficients are calculated using finite difference approximations with 1% changes in the reaction rates :cite:p:`kholodenko1997quantification`.

.. image:: https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/sensitivity_PcFos.png

Control coefficients for integrated pc-Fos are shown by bars (blue, EGF; red, HRG). Numbers above bars indicate the reaction indices, and error bars correspond to simulation standard deviation.


.. note::
    If you want to reuse a result from the previous computation and don't want to calculate sensitivity coefficients again, set ``options['overwrite']`` to ``False``. 