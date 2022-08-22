import numpy as np

from biomass.estimation import convert_scale, initialize_search_param

from .name2idx import C, V
from .ode import initial_values, param_values


class SearchParam(object):
    """Specify model parameters and/or initial values to optimize."""

    def __init__(self):
        # parameters
        self.idx_params = [
            C.V1,
            C.Km1,
            C.V5,
            C.Km5,
            C.V10,
            C.Km10,
            C.n10,
            C.p11,
            C.p12,
            C.p13,
            C.V14,
            C.Km14,
            C.V15,
            C.Km15,
            C.KimDUSP,
            C.KexDUSP,
            C.V20,
            C.Km20,
            C.V21,
            C.Km21,
            C.V24,
            C.Km24,
            C.V25,
            C.Km25,
            C.KimRSK,
            C.KexRSK,
            C.V27,
            C.Km27,
            C.V28,
            C.Km28,
            C.V29,
            C.Km29,
            C.V30,
            C.Km30,
            C.V31,
            C.Km31,
            C.n31,
            C.p32,
            C.p33,
            C.p34,
            C.V35,
            C.Km35,
            C.V36,
            C.Km36,
            C.V37,
            C.Km37,
            C.KimFOS,
            C.KexFOS,
            C.V42,
            C.Km42,
            C.V43,
            C.Km43,
            C.V44,
            C.Km44,
            C.p47,
            C.m47,
            C.p48,
            C.p49,
            C.m49,
            C.p50,
            C.p51,
            C.m51,
            C.V57,
            C.Km57,
            C.n57,
            C.p58,
            C.p59,
            C.p60,
            C.p61,
            C.KimF,
            C.KexF,
            C.p63,
            C.KF31,
            C.nF31,
            C.a,
        ]

        # initial values
        self.idx_initials = [
            # V.(specie)
        ]

    def get_region(self):
        x = param_values()
        y0 = initial_values()

        search_param = initialize_search_param(
            parameters=C.NAMES,
            species=V.NAMES,
            param_values=x,
            initial_values=y0,
            estimated_params=self.idx_params,
            estimated_initials=self.idx_initials,
        )

        search_rgn = np.zeros((2, len(x) + len(y0)))
        # Default: 0.1 ~ 10
        for i, j in enumerate(self.idx_params):
            search_rgn[0, j] = search_param[i] * 0.1  # lower bound
            search_rgn[1, j] = search_param[i] * 10.0  # upper bound
        # Default: 0.5 ~ 2
        for i, j in enumerate(self.idx_initials):
            search_rgn[0, j + len(x)] = search_param[i + len(self.idx_params)] * 0.5  # lower bound
            search_rgn[1, j + len(x)] = search_param[i + len(self.idx_params)] * 2.0  # upper bound

        # search_rgn[:,C.parameter] = [lower_bound, upper_bound]
        # search_rgn[:,V.specie+len(x)] = [lower_bound, upper_bound]

        search_rgn[:, C.V1] = [7.33e-2, 6.60e-01]
        search_rgn[:, C.Km1] = [1.83e2, 8.50e2]
        search_rgn[:, C.V5] = [6.48e-3, 7.20e1]
        search_rgn[:, C.Km5] = [6.00e-1, 1.60e04]
        search_rgn[:, C.V10] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km10] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.n10] = [1.00, 4.00]
        search_rgn[:, C.p11] = [8.30e-13, 1.44e-2]
        search_rgn[:, C.p12] = [8.00e-8, 5.17e-2]
        search_rgn[:, C.p13] = [1.38e-7, 4.84e-1]
        search_rgn[:, C.V14] = [4.77e-3, 4.77e1]
        search_rgn[:, C.Km14] = [2.00e2, 2.00e6]
        search_rgn[:, C.V15] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km15] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.KimDUSP] = [2.20e-4, 5.50e-1]
        search_rgn[:, C.KexDUSP] = [2.60e-4, 6.50e-1]
        search_rgn[:, C.V20] = [4.77e-3, 4.77e1]
        search_rgn[:, C.Km20] = [2.00e2, 2.00e6]
        search_rgn[:, C.V21] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km21] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.V24] = [4.77e-2, 4.77e0]
        search_rgn[:, C.Km24] = [2.00e3, 2.00e5]
        search_rgn[:, C.V25] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km25] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.KimRSK] = [2.20e-4, 5.50e-1]
        search_rgn[:, C.KexRSK] = [2.60e-4, 6.50e-1]
        search_rgn[:, C.V27] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km27] = [1.00e2, 1.00e4]
        search_rgn[:, C.V28] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km28] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.V29] = [4.77e-2, 4.77e0]
        search_rgn[:, C.Km29] = [2.93e3, 2.93e5]
        search_rgn[:, C.V30] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km30] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.V31] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km31] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.n31] = [1.00, 4.00]
        search_rgn[:, C.p32] = [8.30e-13, 1.44e-2]
        search_rgn[:, C.p33] = [8.00e-8, 5.17e-2]
        search_rgn[:, C.p34] = [1.38e-7, 4.84e-1]
        search_rgn[:, C.V35] = [4.77e-3, 4.77e1]
        search_rgn[:, C.Km35] = [2.00e2, 2.00e6]
        search_rgn[:, C.V36] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km36] = [1.00e2, 1.00e4]
        search_rgn[:, C.V37] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km37] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.KimFOS] = [2.20e-4, 5.50e-1]
        search_rgn[:, C.KexFOS] = [2.60e-4, 6.50e-1]
        search_rgn[:, C.V42] = [4.77e-3, 4.77e1]
        search_rgn[:, C.Km42] = [2.00e2, 2.00e6]
        search_rgn[:, C.V43] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km43] = [1.00e2, 1.00e4]
        search_rgn[:, C.V44] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km44] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.p47] = [1.45e-4, 1.45e0]
        search_rgn[:, C.m47] = [6.00e-3, 6.00e1]
        search_rgn[:, C.p48] = [2.70e-3, 2.70e1]
        search_rgn[:, C.p49] = [5.00e-5, 5.00e-1]
        search_rgn[:, C.m49] = [5.00e-3, 5.00e1]
        search_rgn[:, C.p50] = [3.00e-3, 3.00e1]
        search_rgn[:, C.p51] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.m51] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.V57] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.Km57] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.n57] = [1.00, 4.00]
        search_rgn[:, C.p58] = [8.30e-13, 1.44e-2]
        search_rgn[:, C.p59] = [8.00e-8, 5.17e-2]
        search_rgn[:, C.p60] = [1.38e-7, 4.84e-1]
        search_rgn[:, C.p61] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.KimF] = [2.20e-4, 5.50e-1]
        search_rgn[:, C.KexF] = [2.60e-4, 6.50e-1]
        search_rgn[:, C.p63] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.KF31] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.nF31] = [1.00, 4.00]
        search_rgn[:, C.a] = [1.00e2, 5.00e2]

        search_rgn = convert_scale(
            region=search_rgn,
            parameters=C.NAMES,
            species=V.NAMES,
            estimated_params=self.idx_params,
            estimated_initials=self.idx_initials,
        )

        return search_rgn

    def update(self, indiv):
        x = param_values()
        y0 = initial_values()

        for i, j in enumerate(self.idx_params):
            x[j] = indiv[i]
        for i, j in enumerate(self.idx_initials):
            y0[j] = indiv[i + len(self.idx_params)]

        # constraints ----------------------------------------------------------
        x[C.V6] = x[C.V5]
        x[C.Km6] = x[C.Km5]
        x[C.KimpDUSP] = x[C.KimDUSP]
        x[C.KexpDUSP] = x[C.KexDUSP]
        x[C.KimpcFOS] = x[C.KimFOS]
        x[C.KexpcFOS] = x[C.KexFOS]
        x[C.p52] = x[C.p47]
        x[C.m52] = x[C.m47]
        x[C.p53] = x[C.p48]
        x[C.p54] = x[C.p49]
        x[C.m54] = x[C.m49]
        x[C.p55] = x[C.p50]
        x[C.p56] = x[C.p51]
        x[C.m56] = x[C.m51]
        # ----------------------------------------------------------------------

        return x, y0
