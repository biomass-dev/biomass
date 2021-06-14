Usage
=====

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

.. code-block:: python

    >>> from biomass import Model
    >>> from biomass.models import Nakakuki_Cell_2010

    >>> model = Model(Nakakuki_Cell_2010.__package__).create(show_info=True)

    Nakakuki_Cell_2010 information
    ------------------------------
    36 species
    115 parameters, of which 75 to be estimated



Parameter estimation
--------------------

Parameters are adjusted to minimize the distance between model simulation and experimental data.

.. code-block:: python

    from biomass import optimize

    optimize(
        model=model, start=1, options={
            "popsize": 3,
            "max_generation": 100,
            "allowable_error": 0.5,
            "local_search_method": "DE",
            "maxiter": 50,
        }
    )

The temporary result will be saved in ``out/n/`` after each iteration.

Progress list: ``out/n/optimization.log``::

    Generation1: Best Fitness = 5.864228e+00
    Generation2: Best Fitness = 5.864228e+00
    Generation3: Best Fitness = 4.488934e+00
    Generation4: Best Fitness = 3.793744e+00
    Generation5: Best Fitness = 3.652047e+00
    Generation6: Best Fitness = 3.652047e+00
    Generation7: Best Fitness = 3.652047e+00
    Generation8: Best Fitness = 3.452999e+00
    Generation9: Best Fitness = 3.180878e+00
    Generation10: Best Fitness = 1.392501e+00
    Generation11: Best Fitness = 1.392501e+00
    Generation12: Best Fitness = 1.392501e+00
    Generation13: Best Fitness = 1.392501e+00
    Generation14: Best Fitness = 7.018051e-01
    Generation15: Best Fitness = 7.018051e-01
    Generation16: Best Fitness = 7.018051e-01
    Generation17: Best Fitness = 7.018051e-01
    Generation18: Best Fitness = 7.018051e-01
    Generation19: Best Fitness = 6.862063e-01
    Generation20: Best Fitness = 6.862063e-01

* If you want to continue from where you stopped in the last parameter search,

.. code-block:: python

    from biomass import optimize_continue

    optimize_continue(
        model=model, start=1, options={
            "popsize": 3,
            "max_generation": 200,
            "allowable_error": 0.5,
            "local_search_method": "DE",
            "maxiter": 50,
        }
    )

* If you want to search multiple parameter sets (e.g., from 1 to 10) simultaneously,

.. code-block:: python

    from biomass import optimize

    optimize(
        model=model, start=1, end=10, options={
            "popsize": 5,
            "max_generation": 100,
            "allowable_error": 0.5,
            "local_search_method": "DE",
            "maxiter": 50,
        }
    )

* Data Export and Visualization

.. code-block:: python

    from biomass.result import OptimizationResults

    res = OptimizationResults(model)
    # Export estimated parameters in CSV format
    res.to_csv()
    # Visualize estimated parameter sets
    res.savefig(figsize=(16,5), boxplot_kws={"orient": "v"})

.. image:: https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/estimated_parameter_sets.png

Visualization of simulation results
-----------------------------------

.. code-block:: python

    from biomass import run_simulation

    run_simulation(model, viz_type='average', show_all=False, stdev=True)

.. image:: https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/simulation_average.png

Points (blue diamonds, EGF; red squares, HRG) denote experimental data, solid lines denote simulations.

Sensitivity analysis
--------------------

.. code-block:: python

    from biomass import run_analysis

    run_analysis(model, target='reaction', metric='integral', style='barplot')

The single parameter sensitivity of each reaction is defined by

.. math:: C^{M}_{i} = d \ln{M} / d \ln{v_{i}}

where v\ :sub:`i`\ is the i\ :sup:`th`\ reaction rate, v is reaction vector v = (v\ :sub:`1`\, v\ :sub:`2`\, ...) and M is a signaling metric, e.g., time-integrated response, duration.
Sensitivity coefficients were calculated using finite difference approximations with 1% changes in the reaction rates.

.. image:: https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/sensitivity_PcFos.png

Control coefficients for integrated pc-Fos are shown by bars (blue, EGF; red, HRG). Numbers above bars indicate the reaction indices, and error bars correspond to simulation standard deviation.

