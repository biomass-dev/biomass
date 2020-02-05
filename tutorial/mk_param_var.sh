# script to extract parameters from parameter constant
grep '\[C' ../biomass/model/param_const.py | \
sed 's/x\[C\.//g' | \
sed 's/\].*//g' | \
awk '!seen[$0]++' | \
sed -e 's/[[:space:]]\+//g' | \
sed -e "s/\(.*\)/'\1'/" | \
sed 's/$/,/' |
sed -e 's/^/\t/g' > param_var_mid

echo 'param_names = [\' > first_line

echo "    ##
    'len_f_params'\
]

for idx,name in enumerate(param_names):
    exec('%s=%d'%(name,idx))" > last_line


cat first_line param_var_mid last_line > ../biomass/model/name2idx/parameters.py

rm first_line
rm param_var_mid
rm last_line

echo "done processing parameters"

# script to extract variables from variables constant
grep '\[V' ../biomass/model/differential_equation.py | \
sed 's/dydt\[V\.//g' | \
sed 's/\].*//g' | \
awk '!seen[$0]++' | \
sed -e 's/[[:space:]]\+//g' | \
sed -e "s/\(.*\)/'\1'/" | \
sed 's/$/,/' |
sed -e 's/^/\t/g' > var_var_mid

echo 'param_names = [\' > first_line

echo "    ##
    'len_f_params'\
]

for idx,name in enumerate(param_names):
    exec('%s=%d'%(name,idx))" > last_line


cat first_line var_var_mid last_line > ../biomass/model/name2idx/variables.py

rm first_line
rm var_var_mid
rm last_line

echo "done processing variables"