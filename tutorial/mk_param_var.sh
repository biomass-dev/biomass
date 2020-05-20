#!/bin/sh

# script to extract parameters from set_model.py
# extract lines between line containing def f_params and def initial_values
sed -n '/def f_params/,/def initial_values/p' ../biomass/model/set_model.py | \
grep '\[C' | \
sed 's/x\[C\.//g' | \
sed 's/\].*//g' | \
awk '!seen[$0]++' | \
sed -e 's/[[:space:]]\+//g' | \
sed -e "s/\(.*\)/'\1'/" | \
sed 's/$/,/' |
sed -e 's/^/    /g' |
sed -e 's/    //g' > param_var_mid

printf 'param_names = [\n' > first_line

printf "]	

for idx, name in enumerate(param_names):
    exec(
        '{} = {:d}'.format(
            name, idx
        )
    )

len_f_params = len(param_names)" > last_line


cat first_line param_var_mid last_line > ../biomass/model/name2idx/parameters.py

rm first_line
rm param_var_mid
rm last_line


echo "done processing parameters"

# script to extract variables from differential_equation.py
grep 'dydt\[V' ../biomass/model/set_model.py | \
sed 's/dydt\[V\.//g' | \
sed 's/\].*//g' | \
awk '!seen[$0]++' | \
sed -e 's/[[:space:]]\+//g' | \
sed -e "s/\(.*\)/'\1'/" | \
sed 's/$/,/' |
sed -e 's/^/   /g' |
sed -e 's/    //g' > var_var_mid

printf 'var_names = [\n' > first_line

printf "]

for idx, name in enumerate(var_names):
    exec(
        '{} = {:d}'.format(
            name, idx
        )
    )

len_f_vars = len(var_names)" > last_line


cat first_line var_var_mid last_line > ../biomass/model/name2idx/variables.py

rm first_line
rm var_var_mid
rm last_line

echo "done processing variables"