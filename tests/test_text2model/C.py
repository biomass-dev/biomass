INDENT = " " * 4


PPMEK = """\

    @staticmethod
    def _get_ppMEK_slope(t, ligand):
        assert ligand in ['EGF', 'HRG']
        timepoints = [0, 300, 600, 900, 1800, 2700, 3600, 5400]
        ppMEK_data = {
            'EGF': [0.000, 0.773, 0.439, 0.252, 0.130, 0.087, 0.080, 0.066],
            'HRG': [0.000, 0.865, 1.000, 0.837, 0.884, 0.920, 0.875, 0.789],
        }
        assert len(timepoints) == len(ppMEK_data[ligand])
        slope = [
            (ppMEK_data[ligand][i + 1] - activity) / (timepoints[i + 1] - timepoint)
            for i, (timepoint, activity) in enumerate(zip(timepoints, ppMEK_data[ligand]))
            if i + 1 < len(timepoints)
        ]
        for i, timepoint in enumerate(timepoints):
            if timepoint <= t <= timepoints[i + 1]:
                return slope[i]
        assert False

"""


CONDITIONS = """\

        if x[C.Ligand] == 1:  # EGF=10nM
            dydt[V.ppMEKc] = self._get_ppMEK_slope(t, 'EGF')
        elif x[C.Ligand] == 2:  # HRG=10nM
            dydt[V.ppMEKc] = self._get_ppMEK_slope(t, 'HRG')
        else:  # Default: No ligand input
            dydt[V.ppMEKc] = 0.0

"""


BOUNDS = """\

        search_rgn[:, C.V1] = [7.33e-2, 6.60e-01]
        search_rgn[:, C.K1] = [1.83e2, 8.50e2]
        search_rgn[:, C.V5] = [6.48e-3, 7.20e1]
        search_rgn[:, C.K5] = [6.00e-1, 1.60e04]
        search_rgn[:, C.V10] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K10] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.n10] = [1.00, 4.00]
        search_rgn[:, C.kf11] = [8.30e-13, 1.44e-2]
        search_rgn[:, C.kf12] = [8.00e-8, 5.17e-2]
        search_rgn[:, C.kf13] = [1.38e-7, 4.84e-1]
        search_rgn[:, C.V14] = [4.77e-3, 4.77e1]
        search_rgn[:, C.K14] = [2.00e2, 2.00e6]
        search_rgn[:, C.V15] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K15] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.kf18] = [2.20e-4, 5.50e-1]
        search_rgn[:, C.kr18] = [2.60e-4, 6.50e-1]
        search_rgn[:, C.V20] = [4.77e-3, 4.77e1]
        search_rgn[:, C.K20] = [2.00e2, 2.00e6]
        search_rgn[:, C.V21] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K21] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.V24] = [4.77e-2, 4.77e0]
        search_rgn[:, C.K24] = [2.00e3, 2.00e5]
        search_rgn[:, C.V25] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K25] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.kf26] = [2.20e-4, 5.50e-1]
        search_rgn[:, C.kr26] = [2.60e-4, 6.50e-1]
        search_rgn[:, C.V27] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K27] = [1.00e2, 1.00e4]
        search_rgn[:, C.V28] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K28] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.V29] = [4.77e-2, 4.77e0]
        search_rgn[:, C.K29] = [2.93e3, 2.93e5]
        search_rgn[:, C.V30] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K30] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.V31] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K31] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.n31] = [1.00, 4.00]
        search_rgn[:, C.kf32] = [8.30e-13, 1.44e-2]
        search_rgn[:, C.kf33] = [8.00e-8, 5.17e-2]
        search_rgn[:, C.kf34] = [1.38e-7, 4.84e-1]
        search_rgn[:, C.V35] = [4.77e-3, 4.77e1]
        search_rgn[:, C.K35] = [2.00e2, 2.00e6]
        search_rgn[:, C.V36] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K36] = [1.00e2, 1.00e4]
        search_rgn[:, C.V37] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K37] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.kf40] = [2.20e-4, 5.50e-1]
        search_rgn[:, C.kr40] = [2.60e-4, 6.50e-1]
        search_rgn[:, C.V42] = [4.77e-3, 4.77e1]
        search_rgn[:, C.K42] = [2.00e2, 2.00e6]
        search_rgn[:, C.V43] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K43] = [1.00e2, 1.00e4]
        search_rgn[:, C.V44] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K44] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.kf47] = [1.45e-4, 1.45e0]
        search_rgn[:, C.kr47] = [6.00e-3, 6.00e1]
        search_rgn[:, C.kf48] = [2.70e-3, 2.70e1]
        search_rgn[:, C.kf49] = [5.00e-5, 5.00e-1]
        search_rgn[:, C.kr49] = [5.00e-3, 5.00e1]
        search_rgn[:, C.kf50] = [3.00e-3, 3.00e1]
        search_rgn[:, C.kf51] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.kr51] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.V57] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.K57] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.n57] = [1.00, 4.00]
        search_rgn[:, C.kf58] = [8.30e-13, 1.44e-2]
        search_rgn[:, C.kf59] = [8.00e-8, 5.17e-2]
        search_rgn[:, C.kf60] = [1.38e-7, 4.84e-1]
        search_rgn[:, C.kf61] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.kf62] = [2.20e-4, 5.50e-1]
        search_rgn[:, C.kr62] = [2.60e-4, 6.50e-1]
        search_rgn[:, C.kf63] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.KF31] = [np.exp(-10), np.exp(10)]
        search_rgn[:, C.nF31] = [1.00, 4.00]
        search_rgn[:, C.a] = [1.00e2, 5.00e2]

"""

NORMALIZATION = """\
        for observable in self.obs_names:
            self.normalization[observable] = {"timepoint": None, "condition": []}
"""

EXPERIMENTAL_DATA = """\
    def set_data(self):
        self.experiments[self.obs_names.index("Phosphorylated_MEKc")] = {
            "EGF": [0.000, 0.773, 0.439, 0.252, 0.130, 0.087, 0.080, 0.066],
            "HRG": [0.000, 0.865, 1.000, 0.837, 0.884, 0.920, 0.875, 0.789],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_MEKc")] = {
            "EGF": [
                sd / np.sqrt(3) for sd in [0.000, 0.030, 0.048, 0.009, 0.009, 0.017, 0.012, 0.008]
            ],
            "HRG": [
                sd / np.sqrt(3) for sd in [0.000, 0.041, 0.000, 0.051, 0.058, 0.097, 0.157, 0.136]
            ],
        }

        self.experiments[self.obs_names.index("Phosphorylated_ERKc")] = {
            "EGF": [0.000, 0.867, 0.799, 0.494, 0.313, 0.266, 0.200, 0.194],
            "HRG": [0.000, 0.848, 1.000, 0.971, 0.950, 0.812, 0.747, 0.595],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_ERKc")] = {
            "EGF": [
                sd / np.sqrt(3) for sd in [0.000, 0.137, 0.188, 0.126, 0.096, 0.087, 0.056, 0.012]
            ],
            "HRG": [
                sd / np.sqrt(3) for sd in [0.000, 0.120, 0.000, 0.037, 0.088, 0.019, 0.093, 0.075]
            ],
        }

        self.experiments[self.obs_names.index("Phosphorylated_RSKw")] = {
            "EGF": [0, 0.814, 0.812, 0.450, 0.151, 0.059, 0.038, 0.030],
            "HRG": [0, 0.953, 1.000, 0.844, 0.935, 0.868, 0.779, 0.558],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_RSKw")] = {
            "EGF": [
                sd / np.sqrt(3) for sd in [0, 0.064, 0.194, 0.030, 0.027, 0.031, 0.043, 0.051]
            ],
            "HRG": [
                sd / np.sqrt(3) for sd in [0, 0.230, 0.118, 0.058, 0.041, 0.076, 0.090, 0.077]
            ],
        }

        self.experiments[self.obs_names.index("Phosphorylated_cFos")] = {
            "EGF": [0, 0.060, 0.109, 0.083, 0.068, 0.049, 0.027, 0.017],
            "HRG": [0, 0.145, 0.177, 0.158, 0.598, 1.000, 0.852, 0.431],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_cFos")] = {
            "EGF": [
                sd / np.sqrt(3) for sd in [0, 0.003, 0.021, 0.013, 0.016, 0.007, 0.003, 0.002]
            ],
            "HRG": [
                sd / np.sqrt(3) for sd in [0, 0.010, 0.013, 0.001, 0.014, 0.000, 0.077, 0.047]
            ],
        }

        # ----------------------------------------------------------------------

        self.experiments[self.obs_names.index("Phosphorylated_CREBw")] = {
            "EGF": [0, 0.446, 0.030, 0.000, 0.000],
            "HRG": [0, 1.000, 0.668, 0.460, 0.340],
        }
        self.error_bars[self.obs_names.index("Phosphorylated_CREBw")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]],
            "HRG": [sd / np.sqrt(3) for sd in [0, 0.0, 0.0, 0.0, 0.0]],
        }
        # ----------------------------------------------------------------------

        self.experiments[self.obs_names.index("cfos_mRNA")] = {
            "EGF": [0, 0.181, 0.476, 0.518, 0.174, 0.026, 0.000],
            "HRG": [0, 0.353, 0.861, 1.000, 0.637, 0.300, 0.059],
        }
        self.error_bars[self.obs_names.index("cfos_mRNA")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0.017, 0.004, 0.044, 0.004, 0.023, 0.007, 0.008]],
            "HRG": [sd / np.sqrt(3) for sd in [0.017, 0.006, 0.065, 0.044, 0.087, 0.023, 0.001]],
        }
        # ----------------------------------------------------------------------

        self.experiments[self.obs_names.index("cFos_Protein")] = {
            "EGF": [0, 0.078, 0.216, 0.240, 0.320, 0.235],
            "HRG": [0, 0.089, 0.552, 0.861, 1.000, 0.698],
        }
        self.error_bars[self.obs_names.index("cFos_Protein")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0, 0.036, 0.028, 0.056, 0.071, 0.048]],
            "HRG": [sd / np.sqrt(3) for sd in [0, 0.021, 0.042, 0.063, 0.000, 0.047]],
        }

        self.experiments[self.obs_names.index("dusp_mRNA")] = {
            "EGF": [0.000, 0.177, 0.331, 0.214, 0.177, 0.231],
            "HRG": [0.000, 0.221, 0.750, 1.000, 0.960, 0.934],
        }
        self.error_bars[self.obs_names.index("dusp_mRNA")] = {
            "EGF": [sd / np.sqrt(3) for sd in [0.033, 0.060, 0.061, 0.032, 0.068, 0.050]],
            "HRG": [sd / np.sqrt(3) for sd in [0.027, 0.059, 0.094, 0.124, 0.113, 0.108]],
        }

    @staticmethod
    def get_timepoint(obs_name) -> List[int]:
        if obs_name in [
            "Phosphorylated_MEKc",
            "Phosphorylated_ERKc",
            "Phosphorylated_RSKw",
            "Phosphorylated_cFos",
        ]:
            return [0, 300, 600, 900, 1800, 2700, 3600, 5400]  # (Unit: sec.)
        elif obs_name == "Phosphorylated_CREBw":
            return [0, 600, 1800, 3600, 5400]
        elif obs_name == "cfos_mRNA":
            return [0, 600, 1200, 1800, 2700, 3600, 5400]
        elif obs_name in ["cFos_Protein", "dusp_mRNA"]:
            return [0, 900, 1800, 2700, 3600, 5400]
        assert False
"""
