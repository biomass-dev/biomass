## Parameter Estimation of ODE Models
The temporary result will be saved in ```out/n/``` after each iteration.
```python
from biomass import optimize

optimize(Nakakuki_Cell_2010, n)
```
Progress list: ```out/n/optimization.log```
```
Generation1: Best Fitness = 1.726069e+00
Generation2: Best Fitness = 1.726069e+00
Generation3: Best Fitness = 1.726069e+00
Generation4: Best Fitness = 1.645414e+00
Generation5: Best Fitness = 1.645414e+00
Generation6: Best Fitness = 1.645414e+00
Generation7: Best Fitness = 1.645414e+00
Generation8: Best Fitness = 1.645414e+00
Generation9: Best Fitness = 1.645414e+00
Generation10: Best Fitness = 1.645414e+00
Generation11: Best Fitness = 1.645414e+00
Generation12: Best Fitness = 1.645414e+00
Generation13: Best Fitness = 1.645414e+00
Generation14: Best Fitness = 1.645414e+00
Generation15: Best Fitness = 1.645414e+00
Generation16: Best Fitness = 1.249036e+00
Generation17: Best Fitness = 1.171606e+00
Generation18: Best Fitness = 1.171606e+00
Generation19: Best Fitness = 1.171606e+00
Generation20: Best Fitness = 1.171606e+00
```

- If you want to continue from where you stopped in the last parameter search,
```python
from biomass import optimize_continue

optimize_continue(Nakakuki_Cell_2010, n)
```
- If you want to search multiple parameter sets (from *n1* to *n2*) simultaneously,
```python
from biomass import optimize

optimize(Nakakuki_Cell_2010, n1, n2)
```
