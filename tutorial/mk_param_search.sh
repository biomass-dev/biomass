#!/bin/sh

# script to extract the parameter index for optimization
grep '\[C' ../biomass/model/param_const.py | \
sed 's/x\[C\.//g' | \
sed 's/\].*//g' | \
awk '!seen[$0]++' | \
sed -e 's/[[:space:]]\+//g' | \
sed 's/$/,/' |
sed -e 's/^/C./g' |
sed -e 's/^/		/g' |
sed -e 's/    //g' > param_search_mid

# print the header 

echo "import numpy as np

from biomass.model.name2idx import parameters as C
from biomass.model.name2idx import variables as V
from biomass.model.param_const import f_params
from biomass.model.initial_condition import initial_values


def search_parameter_index():

	# Write param index
    search_idx_const = [" > first_line



# last line bugged, therefore default will be used

cat first_line param_search_mid last_line_paramsearch > ../biomass/param_estim/search_parameter.py

rm first_line
rm param_search_mid


echo "done processing search param"
