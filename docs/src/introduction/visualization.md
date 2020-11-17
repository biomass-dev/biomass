## Visualization of Simulation Results
In BioMASS, you can visualize your model using an optimized parameter set, or simulate your model with user-defined values for parameters and initial amounts.

## Customizing your visualization
In the *mapk_cascade* model directory, open the file `viz.py`. It contains all necessary and customizable information to adapt your output plots.
By default, you can observe multiple plots each showing a single observable, or plot several observables in a single plot (multiplot). Note that you cannot generate several multiplots in one simulation. Instead, you need to simulate your model sequentially as many times as your desired number of multiplots (commenting out all other multiplots).

Without changing anything in the `viz.py` file, BioMASS will plot the observables which you defined in the `observables.py` file. If you would like to customize the y-axis label for your species name, you can do this within the `get_timecourse_options(self)` function starting from line 112. In the *mapk_cascade* model, the plot label for 'biphosphorylated_MAPK' has been changed to 'ERK-PP', and, similarly, 'unphosphorylated_MAPK' to 'ERK'. 

## Running a simulation and visualize the results
We need to simulate and visualize the model from the biomass folder. 

In your terminal, start Python from the biomass folder:
`$ python3
`
<br>
Next, load the required modules: <br>
`
$ >>> from biomass.models import mapk_cascade
`
<br>
`
$ >>> from biomass import run_simulation
`

Let's visualize the results from our previous parameter optimization. Since we only calculated one parameter set, we can only choose the "best" option to do so. If you optimize more than one parameter set, you can also show the average of your parameter sets. <br>
`
$ >>> run_simulation(mapk_cascade, viz_type='best', show_all=False, stdev=True)
`

The plots have been placed in the automatically generated directory "best", and can be found in:
<br>
`biomass->biomass-> models-> mapk_cascade -> figure -> simulation -> best`
<br>

Note: If you run Python from the terminal and you update your `viz.py` file, you need to quit Python, restart it, and run the above commands again to put the changes into effect. 
Alternatively, you can save the three python commands into one file and place it into the biomass folder. You can then run your simulation with:

`
$ python3 <simulation_file_name.py>
`

Replace `<simulation_file_name.py>` with the name of the file you previously created.


