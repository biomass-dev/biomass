from .name2idx import C, V
from .reaction_network import ReactionNetwork


class DifferentialEquation(ReactionNetwork):
    def __init__(self, perturbation):
        super(DifferentialEquation, self).__init__()
        self.perturbation = perturbation

    @staticmethod
    def _get_ppMEK_slope(t, ligand) -> float:
        assert ligand in ["EGF", "HRG"]
        timepoints = [0, 300, 600, 900, 1800, 2700, 3600, 5400]
        ppMEK_data = {
            "EGF": [0.000, 0.773, 0.439, 0.252, 0.130, 0.087, 0.080, 0.066],
            "HRG": [0.000, 0.865, 1.000, 0.837, 0.884, 0.920, 0.875, 0.789],
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

    # Refined Model
    def diffeq(self, t, y, *x):

        v = self.flux(t, y, x)

        if self.perturbation:
            for i, dv in self.perturbation.items():
                v[i] = v[i] * dv

        dydt = [0] * V.NUM

        if x[C.Ligand] == x[C.EGF]:  # EGF=10nM
            dydt[V.ppMEKc] = self._get_ppMEK_slope(t, "EGF")
        elif x[C.Ligand] == x[C.HRG]:  # HRG=10nM
            dydt[V.ppMEKc] = self._get_ppMEK_slope(t, "HRG")
        else:  # Default: No ligand input
            dydt[V.ppMEKc] = 0.0

        dydt[V.CREBn] = -v[27] + v[28]
        dydt[V.pCREBn] = v[27] - v[28]
        dydt[V.ERKc] = -v[1] + v[3] - v[7]
        dydt[V.ERKn] = v[5] + v[7] * (x[C.Vc] / x[C.Vn]) + v[50] - v[51] + v[55] - v[56]
        dydt[V.pERKc] = v[1] - v[2] - v[3] + v[4] - v[8]
        dydt[V.pERKn] = -v[5] + v[6] + v[8] * (x[C.Vc] / x[C.Vn]) + v[48] - v[49] + v[53] - v[54]
        dydt[V.ppERKc] = v[2] - v[4] - v[9]
        dydt[V.ppERKn] = -v[6] + v[9] * (x[C.Vc] / x[C.Vn]) - v[47] - v[52]
        dydt[V.Elk1n] = -v[29] + v[30]
        dydt[V.pElk1n] = v[29] - v[30]
        dydt[V.cFOSc] = v[34] - v[35] - v[36] + v[37] - v[38] - v[40]
        dydt[V.cFOSn] = v[40] * (x[C.Vc] / x[C.Vn]) - v[42] - v[43] + v[44] - v[45]
        dydt[V.pcFOSc] = v[35] + v[36] - v[37] - v[39] - v[41]
        dydt[V.pcFOSn] = v[41] * (x[C.Vc] / x[C.Vn]) + v[42] + v[43] - v[44] - v[46]
        dydt[V.DUSPc] = v[13] - v[14] + v[15] - v[16] - v[18]
        dydt[V.DUSPn] = (
            v[18] * (x[C.Vc] / x[C.Vn])
            - v[20]
            + v[21]
            - v[22]
            - v[47]
            + v[48]
            - v[49]
            + v[50]
            - v[51]
        )
        dydt[V.pDUSPc] = v[14] - v[15] - v[17] - v[19]
        dydt[V.pDUSPn] = (
            v[19] * (x[C.Vc] / x[C.Vn])
            + v[20]
            - v[21]
            - v[23]
            - v[52]
            + v[53]
            - v[54]
            + v[55]
            - v[56]
        )
        dydt[V.DUSPn_ERKn] = v[51]
        dydt[V.DUSPn_pERKn] = v[49] - v[50]
        dydt[V.DUSPn_ppERKn] = v[47] - v[48]
        dydt[V.pDUSPn_ERKn] = v[56]
        dydt[V.pDUSPn_pERKn] = v[54] - v[55]
        dydt[V.pDUSPn_ppERKn] = v[52] - v[53]
        dydt[V.RSKc] = -v[24] + v[25]
        dydt[V.pRSKc] = v[24] - v[25] - v[26]
        dydt[V.pRSKn] = v[26] * (x[C.Vc] / x[C.Vn])
        dydt[V.PrecfosmRNAn] = v[31] - v[32]
        dydt[V.PreduspmRNAn] = v[10] - v[11]
        dydt[V.cfosmRNAc] = v[32] * (x[C.Vn] / x[C.Vc]) - v[33]
        dydt[V.duspmRNAc] = v[11] * (x[C.Vn] / x[C.Vc]) - v[12]
        dydt[V.Fc] = v[60] - v[61] - v[62]
        dydt[V.Fn] = v[62] * (x[C.Vc] / x[C.Vn]) - v[63]
        dydt[V.FmRNAc] = v[58] * (x[C.Vn] / x[C.Vc]) - v[59]
        dydt[V.PreFmRNAn] = v[57] - v[58]

        return dydt


def param_values():

    x = [0] * C.NUM

    x[C.V1] = 0.34284837
    x[C.Km1] = 307.0415253
    x[C.V2] = 2.20e-01
    x[C.Km2] = 3.50e02
    x[C.V3] = 7.20e-01
    x[C.Km3] = 1.60e02
    x[C.V4] = 6.48e-01
    x[C.Km4] = 6.00e01
    x[C.V5] = 19.49872346
    x[C.Km5] = 29.94073716
    x[C.V6] = x[C.V5]
    x[C.Km6] = x[C.Km5]
    x[C.KimERK] = 1.20e-02
    x[C.KexERK] = 1.80e-02
    x[C.KimpERK] = 1.20e-02
    x[C.KexpERK] = 1.80e-02
    x[C.KimppERK] = 1.10e-02
    x[C.KexppERK] = 1.30e-02
    x[C.V10] = 29.24109258
    x[C.Km10] = 169.0473748
    x[C.n10] = 3.970849295
    x[C.p11] = 0.000126129
    x[C.p12] = 0.007875765
    x[C.p13] = 0.001245747
    x[C.V14] = 5.636949216
    x[C.Km14] = 34180.48
    x[C.V15] = 2.992346912
    x[C.Km15] = 0.001172165
    x[C.p16] = 2.57e-04
    x[C.p17] = 9.63e-05
    x[C.KimDUSP] = 0.024269764
    x[C.KexDUSP] = 0.070467899
    x[C.KimpDUSP] = x[C.KimDUSP]
    x[C.KexpDUSP] = x[C.KexDUSP]
    x[C.V20] = 0.157678678
    x[C.Km20] = 735598.6967
    x[C.V21] = 0.005648117
    x[C.Km21] = 387.8377182
    x[C.p22] = 2.57e-04
    x[C.p23] = 9.63e-05
    x[C.V24] = 0.550346114
    x[C.Km24] = 29516.06587
    x[C.V25] = 10.09063736
    x[C.Km25] = 0.913939859
    x[C.KimRSK] = 0.025925065
    x[C.KexRSK] = 0.129803956
    x[C.V27] = 19.23118154
    x[C.Km27] = 441.5834425
    x[C.V28] = 6.574759504
    x[C.Km28] = 14.99180922
    x[C.V29] = 0.518529841
    x[C.Km29] = 21312.69109
    x[C.V30] = 13.79479021
    x[C.Km30] = 15.04396629
    x[C.V31] = 0.655214248
    x[C.Km31] = 185.9760682
    x[C.n31] = 1.988003164
    x[C.KF31] = 0.013844393
    x[C.nF31] = 2.800340453
    x[C.p32] = 0.003284434
    x[C.p33] = 0.000601234
    x[C.p34] = 7.65e-05
    x[C.V35] = 8.907637012
    x[C.Km35] = 8562.744184
    x[C.V36] = 0.000597315
    x[C.Km36] = 528.552427
    x[C.V37] = 1.745848179
    x[C.Km37] = 0.070379236
    x[C.p38] = 2.57e-04
    x[C.p39] = 9.63e-05
    x[C.KimFOS] = 0.54528521
    x[C.KexFOS] = 0.133249762
    x[C.KimpcFOS] = x[C.KimFOS]
    x[C.KexpcFOS] = x[C.KexFOS]
    x[C.V42] = 0.909968714
    x[C.Km42] = 3992.061328
    x[C.V43] = 0.076717457
    x[C.Km43] = 1157.116021
    x[C.V44] = 0.078344305
    x[C.Km44] = 0.051168202
    x[C.p45] = 2.57e-04
    x[C.p46] = 9.63e-05
    x[C.p47] = 0.001670815
    x[C.m47] = 15.80783969
    x[C.p48] = 0.686020478
    x[C.p49] = 0.314470502
    x[C.m49] = 2.335459127
    x[C.p50] = 26.59483436
    x[C.p51] = 0.01646825
    x[C.m51] = 9.544308421
    x[C.p52] = x[C.p47]
    x[C.m52] = x[C.m47]
    x[C.p53] = x[C.p48]
    x[C.p54] = x[C.p49]
    x[C.m54] = x[C.m49]
    x[C.p55] = x[C.p50]
    x[C.p56] = x[C.p51]
    x[C.m56] = x[C.m51]
    x[C.V57] = 1.026834758
    x[C.Km57] = 0.637490056
    x[C.n57] = 3.584464176
    x[C.p58] = 0.000270488
    x[C.p59] = 0.001443889
    x[C.p60] = 0.002448164
    x[C.p61] = 3.50e-05
    x[C.KimF] = 0.019898797
    x[C.KexF] = 0.396950616
    x[C.p63] = 4.13e-05

    x[C.a] = 218.6276381
    x[C.Vn] = 0.22
    x[C.Vc] = 0.94

    x[C.EGF] = 0
    x[C.HRG] = 1
    x[C.no_ligand] = 2

    return x


def initial_values():

    y0 = [0] * V.NUM

    y0[V.ERKc] = 9.60e02
    y0[V.RSKc] = 3.53e02
    y0[V.CREBn] = 1.00e03
    y0[V.Elk1n] = 1.51e03

    return y0
