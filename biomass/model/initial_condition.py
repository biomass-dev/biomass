from .name2idx import variables as V

def initial_values():
    y0 = [0]*V.len_f_vars

    y0[V.ERKc] = 9.60e+02
    y0[V.RSKc] = 3.53e+02
    y0[V.CREBn] = 1.00e+03
    y0[V.Elk1n] = 1.51e+03

    return y0