# BioMASS
## Modeling and Analysis of Signaling Systems
<img align="left" src=public/images/logo.png width="300">
Mathematical modeling has the power to simplify and understand the input-output mechanism of biological systems. Although there are many researches devoted on producing models to describe dynamical cellular signaling systems, most of these models are limited and do not cover multiple pathways. Therefore, there is a challenge to combine these models to enable understanding at a larger scale. Nevertheless, larger network means that it gets more difficult to estimate parameters to reproduce dynamic experimental data needed for deeper understanding of a system.<br>
To overcome these problems, we developed BioMASS, a modeling platform tailored to optimizing mathematical models of biological processes. By using BioMASS, users can efficiently optimize kinetic parameters to fit user-defined models to experimental data, while performing analysis on reaction networks to predict critical components affecting cellular output.

## Description
BioMASS is a biological modeling environment tailored to

1. Parameter Estimation of ODE Models
1. Sensitivity Analysis

currently implimented for modeling early transcriptional regulation pathway ([Nakakuki *et al., Cell*, 2010](https://doi.org/10.1016/j.cell.2010.03.054)).

## Dependencies
> - Python 3+
> - numpy, scipy, matplotlib, seaborn

## Usage
#### Parameter Estimation of ODE Models (*n* = 1, 2, 3, · · ·)
The temporary result will be saved in **out/*n*/out.log** after each iteration.
```bash
$ nohup python run_ga.py n &
```
- If you want to continue from where you stopped in the last parameter search,
```bash
$ nohup python run_ga_continue.py n &
```
- If you want to search multiple parameter sets (from *n1* to *n2*) simutaneously,
```bash
$ nohup python run_ga.py n1 n2 &
```

---
#### Visualization of Simulation Results
```bash
$ python run_sim.py
```

|viz_type|Description|
|--------|-----------|
|'average'|The average of simulation results with parameter sets in ```out/```|
|'best'|The best simulation result in ```out/```, simulation with ```best_fit_param```|
|'original'|Simulation with the default parameters and initial values defined in ```biomass/model/```|
|'n(=1,2,...)'|Use the parameter set in ```out/n/```|

- ```viz_type='average',show_all=False,stdev=True```
![simulation_average](public/images/simulation_average.png)

- ```viz_type='best',show_all=True,stdev=False```
![simulation_best](public/images/simulation_best.png)

    Points (blue diamonds, EGF; red squares, HRG) denote experimental data, solid lines denote simulations

---
#### Sensitivity Analysis
The single parameter sensitivity of each reaction is defined by<br>

*s<sub>i</sub>*(*q*(**v**),*v<sub>i</sub>*) = *∂* ln(*q*(**v**)) / *∂* ln(*v<sub>i</sub>*) = *∂*_q_(**v**) / *∂*_v<sub>i</sub>_ · *v<sub>i</sub>* / *q*(**v**)

where *v<sub>i</sub>* is the *i*<sup>th</sup> reaction rate, **v** is reaction vector **v** = (*v<sub>1</sub>*, *v<sub>2</sub>*, ...) and *q*(**v**) is a target function, e.g., time-integrated response, duration. Sensitivity coefficients were calculated using finite difference approximations with 1% changes in the reaction rates.

```bash
$ python run_analysis.py
```
![sensitivity_PcFos](public/images/sensitivity_PcFos.png)

Control coefficients for integrated pc-Fos are shown by bars (blue, EGF; red, HRG). Numbers above bars indicate the reaction indices, and error bars correspond to simulation standard deviation.

## Installation
    $ git clone https://github.com/okadalabipr/biomass.git

## References
- Nakakuki, T. *et al.* Ligand-specific c-Fos expression emerges from the spatiotemporal control of ErbB network dynamics. *Cell* **141**, 884–896 (2010). https://doi.org/10.1016/j.cell.2010.03.054

- Kimura, S., Ono, I., Kita, H. & Kobayashi, S. An extension of UNDX based on guidelines for designing crossover operators: proposition and evaluation of ENDX. *Trans. Soc. Instrum. Control Eng.* **36**, 1162–1171 (2000). https://doi.org/10.9746/sicetr1965.36.1162

- Kimura, S. & Konagaya, A. A Genetic Algorithm with Distance Independent Diversity Control for High Dimensional Function Optimization. *J. Japanese Soc. Artif. Intell.* **18**, 193–202 (2003). https://doi.org/10.1527/tjsai.18.193

- Kimura, S., Nakakuki, T., Kirita, S. & Okada, M. AGLSDC: A Genetic Local Search Suitable for Parallel Computation. *SICE J. Control. Meas. Syst. Integr.* **4**, 105–113 (2012). https://doi.org/10.9746/jcmsi.4.105

- Kholodenko, B. N., Demin, O. V. & Westerhoff, H. V. Control Analysis of Periodic Phenomena in Biological Systems. *J. Phys. Chem. B* **101**, 2070–2081 (1997). https://doi.org/10.1021/jp962336u

- Kholodenko, B. N., Hoek, J. B., Westerhoff, H. V. & Brown, G. C. Quantification of information transfer via cellular signal transduction pathways. *FEBS Lett.* **414**, 430–434 (1997). https://doi.org/10.1016/S0014-5793(97)01018-1

## License
[MIT](/LICENSE)