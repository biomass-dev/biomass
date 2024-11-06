from typing import Dict, List

import numpy as np

from biomass.dynamics.solver import solve_ode

from .name2idx import C, V
from .ode import DifferentialEquation


class Observable(DifferentialEquation):
    """
    Correlating model simulations and experimental measurements.

    Attributes
    ----------
    obs_names : list of strings
        Names of self.obs_names.

    t : range
        Simulation time span.

    conditions : list of strings
        Experimental conditions.

    simulations : numpy.ndarray
        The numpy array to store simulation results.

    normalization : nested dict
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If :obj:`None`, the maximum value will be used for normalization.

        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in ``sim.conditions`` will be used.

    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """

    def __init__(self):
        super(Observable, self).__init__(perturbation={})
        self.obs_names = [
            "IL13_stimulation",
            "RecSurf",
            "IL13-cell",
            "pIL4Ra",
            "pJAK2",
            "SOCS3mRNA",
            "CD274mRNA",
            "SOCS3",
            "pSTAT5",
        ]
        self.t: range = range(0, 120 + 1)
        self.conditions: list = [
            "IL13_0",
            "IL13_4",
            "IL13_20",
            "IL13_80",
        ]
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.conditions), len(self.t))
        )
        self.normalization: dict = {}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

    def simulate(self, x, y0, _perturbation=None):
        if _perturbation is not None:
            self.perturbation = _perturbation

        for i, condition in enumerate(self.conditions):
            if condition == "IL13_0":
                y0[V.IL13] = 0.0
            elif condition == "IL13_4":
                y0[V.IL13] = 4.0
            elif condition == "IL13_20":
                y0[V.IL13] = 20.0
            elif condition == "IL13_80":
                y0[V.IL13] = 80.0

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index("IL13_stimulation"), i] = sol.y[V.IL13]
                self.simulations[self.obs_names.index("RecSurf"), i] = (
                    sol.y[V.Rec] + sol.y[V.IL13_Rec] + sol.y[V.p_IL13_Rec]
                )
                self.simulations[self.obs_names.index("IL13-cell"), i] = (
                    sol.y[V.IL13_Rec]
                    + sol.y[V.p_IL13_Rec]
                    + sol.y[V.p_IL13_Rec_i]
                    + sol.y[V.IL13_DecoyR]
                ) * x[C.scale_IL13_cell]
                self.simulations[self.obs_names.index("pIL4Ra"), i] = (
                    sol.y[V.p_IL13_Rec] + sol.y[V.p_IL13_Rec_i]
                ) * x[C.scale_pIL4Ra]
                self.simulations[self.obs_names.index("pJAK2"), i] = (
                    sol.y[V.pJAK2] * x[C.scale_pJAK2]
                )
                self.simulations[self.obs_names.index("SOCS3mRNA"), i] = (
                    sol.y[V.SOCS3mRNA] * x[C.scale_SOCS3mRNA]
                )
                self.simulations[self.obs_names.index("CD274mRNA"), i] = (
                    sol.y[V.CD274mRNA] * x[C.scale_CD274mRNA]
                )
                self.simulations[self.obs_names.index("SOCS3"), i] = sol.y[V.SOCS3]
                self.simulations[self.obs_names.index("pSTAT5"), i] = sol.y[V.pSTAT5]
        return None

    def set_data(self):
        self.experiments[self.obs_names.index("RecSurf")] = {
            "IL13_0": [
                1.289752779,
                1.57854438,
                1.515091272,
                1.50496416,
                1.617802874,
                1.572371701,
                1.430619477,
            ],
            "IL13_4": [
                1.40469351,
                1.36298049,
                1.333586237,
                1.210536056,
                1.256908767,
                1.255595571,
                1.230653226,
            ],
            "IL13_80": [
                1.321308168,
                1.179069779,
                1.042164123,
                1.056526564,
                1.056091246,
                1.078735295,
                0.929233191,
            ],
        }
        self.experiments[self.obs_names.index("IL13-cell")] = {
            "IL13_0": [0.616, 0.296, 0.322, 0.170, 0.258, 0.232],
            "IL13_4": [
                1.186,
                1.458,
                1.797,
                1.783,
                2.295,
                2.470,
                2.515,
                2.656,
                2.956,
                3.199,
                2.621,
            ],
            "IL13_20": [
                3.975,
                4.130,
                5.016,
                4.933,
                6.025,
                6.352,
                5.879,
                5.076,
                6.211,
                5.947,
                5.095,
            ],
            "IL13_80": [
                9.782,
                9.952,
                10.674,
                11.248,
                12.747,
                12.398,
                10.969,
                8.286,
                11.345,
                11.198,
                10.552,
            ],
        }
        self.experiments[self.obs_names.index("pIL4Ra")] = {
            "IL13_4": [
                0.338107913,
                1.084337349,
                0.539111503,
                0.299787385,
                0.845156744,
                1.175990856,
                1.35788895,
                1.442254159,
                1.123066049,
                1.551873086,
                1.147885063,
                0.995444907,
                1.242317637,
                1.28102483,
                1.194533888,
                1.133331567,
            ],
            "IL13_20": [
                0.185418978,
                0.442239546,
                1.63082953,
                1.926293409,
                1.40013759,
                2.223954642,
                2.258737723,
                1.541743723,
                1.659951486,
                2.283486889,
                1.364989369,
            ],
        }
        self.experiments[self.obs_names.index("pJAK2")] = {
            "IL13_4": [
                0.1871,
                0.5336,
                0.6451,
                0.6394,
                0.2452,
                0.1525,
                0.8347,
                1.4584,
                1.2889,
                1.2754,
                1.4329,
                0.8428,
                1.3672,
                1.7532,
                1.4344,
                1.037,
                1.1503,
                1.3739,
                1.2167,
                1.49,
                0.3656,
                0.6675,
            ],
            "IL13_20": [
                0.2214,
                0.6896,
                1.2783,
                1.6887,
                1.3888,
                2.0209,
                2.1026,
                1.6429,
                2.0556,
                1.7644,
                1.8168,
                1.0997,
                1.4807,
                1.584,
                0.9229,
                1.2044,
            ],
        }
        self.experiments[self.obs_names.index("SOCS3mRNA")] = {
            "IL13_0": [
                1.562357019,
                0.992868163,
                14.83449477,
                11.675177,
                11.62681448,
                5.710282463,
                1.243430459,
                0.659244766,
                14.54733886,
            ],
            "IL13_20": [
                1.399863921,
                15.59376831,
                29.66110989,
                250.7513584,
                208.0216845,
                719.6189149,
                1154.942401,
                504.6313426,
                1017.019576,
            ],
        }
        self.experiments[self.obs_names.index("CD274mRNA")] = {
            "IL13_0": [0.0, -0.105611482, 0.042114502, -0.038136714, 0.009833619],
            "IL13_20": [0.0, -0.1909, 0.6748, 2.305, 1.775],
        }
        self.experiments[self.obs_names.index("SOCS3")] = {
            "IL13_0": [8.33, 6.44, 7.62, 8.69, 13.3, 16.43],
            "IL13_20": [35.1, 23.88, 79.72, 148.68, 259.11, 271.91, 255.38, 204.46],
            "IL13_80": [6.05, 6.23, 87.45, 103.27, 221.74, 262.85, 140.77],
        }
        self.experiments[self.obs_names.index("pSTAT5")] = {
            "IL13_4": [
                1.598406966,
                2.012148639,
                7.644693901,
                30.339995,
                34.49644889,
                54.73958911,
                48.86794422,
                57.76829847,
                68.400004,
                57.5319745,
                84.08635186,
                80.47435164,
                66.61604471,
                57.65815916,
                52.79857246,
                56.17671662,
            ],
            "IL13_20": [
                1.716733095,
                9.772025181,
                32.81159437,
                50.80322648,
                68.99256124,
                107.9066904,
                111.5714137,
                106.7423055,
                104.28,
                125.2740088,
                115.7026455,
                128.2373073,
                123.2374928,
                89.65411083,
                65.83056831,
                108.2978696,
                80.75011175,
            ],
        }

    def get_timepoint(self, obs_name) -> Dict[str, List[int]]:
        if obs_name in self.obs_names:
            if obs_name == "RecSurf":
                return {
                    condition: [0.0, 10.0, 20.0, 30.0, 60.0, 90.0, 120.0]
                    for condition in ["IL13_0", "IL13_4", "IL13_80"]
                }
            elif obs_name == "IL13-cell":
                return {
                    "IL13_0": [0.0, 15.0, 30.0, 60.0, 90.0, 120.0],
                    "IL13_4": [5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0],
                    "IL13_20": [5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0],
                    "IL13_80": [5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0],
                }
            elif obs_name == "pIL4Ra":
                return {
                    "IL13_4": [
                        0.0,
                        4.0,
                        5.0,
                        7.0,
                        10.0,
                        15.0,
                        17.5,
                        20.0,
                        22.5,
                        25.0,
                        30.0,
                        35.0,
                        40.0,
                        45.0,
                        50.0,
                        60.0,
                    ],
                    "IL13_20": [0.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0, 45.0, 60.0],
                }
            elif obs_name == "pJAK2":
                return {
                    "IL13_4": [
                        0.0,
                        2.5,
                        4.0,
                        5.0,
                        7.0,
                        7.5,
                        10.0,
                        12.5,
                        15.0,
                        17.5,
                        20.0,
                        25.0,
                        30.0,
                        35.0,
                        40.0,
                        45.0,
                        50.0,
                        60.0,
                        75.0,
                        90.0,
                        105.0,
                        120.0,
                    ],
                    "IL13_20": [
                        0.0,
                        4.0,
                        5.0,
                        7.0,
                        10.0,
                        15.0,
                        20.0,
                        25.0,
                        30.0,
                        40.0,
                        45.0,
                        50.0,
                        60.0,
                        80.0,
                        100.0,
                        120.0,
                    ],
                }
            elif obs_name == "SOCS3mRNA":
                return {
                    condition: [0.0, 20.0, 30.0, 40.0, 60.0, 80.0, 90.0, 100.0, 120.0]
                    for condition in ["IL13_0", "IL13_20"]
                }
            elif obs_name == "CD274mRNA":
                return {
                    condition: [0.0, 30.0, 60.0, 90.0, 120.0]
                    for condition in ["IL13_0", "IL13_20"]
                }
            elif obs_name == "SOCS3":
                return {
                    "IL13_0": [0.0, 20.0, 40.0, 60.0, 90.0, 120.0],
                    "IL13_20": [0.0, 20.0, 40.0, 60.0, 80.0, 90.0, 100.0, 120.0],
                    "IL13_80": [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0],
                }
            elif obs_name == "pSTAT5":
                return {
                    "IL13_4": [
                        0.0,
                        5.0,
                        10.0,
                        15.0,
                        20.0,
                        25.0,
                        30.0,
                        35.0,
                        40.0,
                        45.0,
                        50.0,
                        60.0,
                        70.0,
                        80.0,
                        100.0,
                        120.0,
                    ],
                    "IL13_20": [
                        0.0,
                        5.0,
                        10.0,
                        15.0,
                        20.0,
                        25.0,
                        30.0,
                        35.0,
                        40.0,
                        45.0,
                        50.0,
                        60.0,
                        70.0,
                        80.0,
                        90.0,
                        100.0,
                        120.0,
                    ],
                }
        assert False
