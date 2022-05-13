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
    
    optimize(model, x_id=1)

The temporary result will be saved in ``out/_tmp/`` after each iteration.

Progress list: ``out/_tmp/optimization.log``::

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
    
.. warning::

    To set optimizer_options["workers"] greater than 1, use :class:`~biomass.estimation.ExternalOptimizer` (see example below).

* If you want to search multiple parameter sets (e.g., from 1 to 10) simultaneously,

1. Prepare ``optimize.py``

.. code-block:: python
    
    import sys
    from biomass import Model
    from biomass.models import Nakakuki_Cell_2010
    
    if __name__ == "__main__":
        args = sys.argv
        model = Model(Nakakuki_Cell_2010.__package__).create()
        optimize(model, x_id=args[1], disp_here=True)

2. Prepare ``optimize_parallel.sh``

.. code-block:: shell
    
    #!/bin/sh
    
    for i in $(seq 1 10); do
        nohup python optimzie.py $i >> progress/$i.log 2>&1 &
    done

3. Run ``optimize_parallel.sh``

.. code-block::
    
    $ mkdir progress
    $ sh optimize_parallel.sh

To kill jobs, run

.. code-block::
    
    $ pgrep -f optimize.py | xargs kill -9

Using external optimizers
^^^^^^^^^^^^^^^^^^^^^^^^^

You can also use external optimization methods to determine model parameters.
Below is an example of using ``scipy.optimize.differential_evolution`` for parameter estimation.

.. code-block:: python

    from scipy.optimize import differential_evolution

    from biomass import Model, run_simulation
    from biomass.estimation import Optimizer
    from biomass.models import Nakakuki_Cell_2010

    model = Model(Nakakuki_Cell_2010.__package__).create()
    param_idx = 1
    optimizer = Optimizer(model, differential_evolution, param_idx)

    def obj_fun(x):
        """Objective function to be minimized."""
        return optimizer.get_obj_val(x)

    res = optimizer.minimize(
        obj_fun,
        [(0, 1) for _ in range(len(model.problem.bounds))],
        strategy="best1bin",
        maxiter=50,
        tol=1e-4,
        mutation=0.1,
        recombination=0.5,
        disp=True,
        polish=False,
        workers=-1,
    )
    
    # Import the solution of the optimization (res.x) and visualize the result.
    param_values = model.problem.gene2val(res.x)
    optimizer.import_solution(param_values)
    run_simulation(model, viz_type=str(param_idx))

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