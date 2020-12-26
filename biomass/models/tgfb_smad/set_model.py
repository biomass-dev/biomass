from .name2idx import C, V


class DifferentialEquation(object):
    def __init__(self, perturbation):
        self.perturbation = perturbation

    def diffeq(self, t, y, *x):

        v = {}

        v[1] = y[V.Rec] * x[C.Rec_act] * y[V.TGFb]
        v[2] = y[V.TGFb_pRec] * x[C.pRec_degind]
        v[3] = x[C.k_on_222] * (y[V.ppS2] ** 3)
        v[4] = 3 * x[C.S_dephosphos] * y[V.ppS2_ppS2_ppS2]
        v[5] = x[C.k_on_333] * y[V.ppS3] ** 3
        v[6] = 3 * x[C.S_dephosphos] * y[V.ppS3_ppS3_ppS3]
        v[7] = (y[V.S4] ** 3) * x[C.k_on_444]
        v[8] = y[V.S4_S4_S4] * x[C.kdiss_SS]
        v[9] = y[V.S2] * x[C.S_phos] * y[V.TGFb_pRec]
        v[10] = x[C.S_dephosphos] * y[V.ppS2]
        v[11] = x[C.S_dephos] * y[V.pS2]
        v[12] = y[V.S3] * x[C.S_phos] * y[V.TGFb_pRec]
        v[13] = x[C.S_dephosphos] * y[V.ppS3]
        v[14] = x[C.S_dephos] * y[V.pS3]
        v[15] = x[C.k_on_223] * (y[V.ppS2] ** 2) * y[V.ppS3]
        v[16] = 2 * x[C.S_dephosphos] * y[V.ppS2_ppS2_ppS3]
        v[17] = x[C.S_dephosphos] * y[V.ppS2_ppS2_ppS3]
        v[18] = y[V.S4] * x[C.k_on_224] * (y[V.ppS2] ** 2)
        v[19] = 2 * x[C.S_dephosphos] * y[V.ppS2_ppS2_S4]
        v[20] = x[C.k_on_233] * y[V.ppS2] * (y[V.ppS3] ** 2)
        v[21] = x[C.S_dephosphos] * y[V.ppS2_ppS3_ppS3]
        v[22] = 2 * x[C.S_dephosphos] * y[V.ppS2_ppS3_ppS3]
        v[23] = y[V.S4] * x[C.k_on_334] * (y[V.ppS3] ** 2)
        v[24] = 2 * x[C.S_dephosphos] * y[V.ppS3_ppS3_S4]
        v[25] = (y[V.S4] ** 2) * x[C.k_on_244] * y[V.ppS2]
        v[26] = x[C.S_dephosphos] * y[V.ppS2_S4_S4]
        v[27] = (y[V.S4] ** 2) * x[C.k_on_344] * y[V.ppS3]
        v[28] = x[C.S_dephosphos] * y[V.ppS3_S4_S4]
        v[29] = y[V.S4] * x[C.k_on_234] * y[V.ppS2] * y[V.ppS3]
        v[30] = x[C.S_dephosphos] * y[V.ppS2_ppS3_S4]
        v[31] = x[C.S_dephosphos] * y[V.ppS2_ppS3_S4]
        v[32] = (
            x[C.gene_turn]
            + x[C.gene_act1] * y[V.ppS2_ppS3_ppS3]
            + x[C.gene_act2] * y[V.ppS2_S4_S4]
            + x[C.gene_act3] * y[V.ppS2_ppS3_S4]
        ) / (
            x[C.gene_inh1] * y[V.ppS2_ppS3_ppS3]
            + x[C.gene_inh2] * y[V.ppS2_S4_S4]
            + x[C.gene_inh3] * y[V.ppS2_ppS3_S4]
            + 1
        )
        v[33] = y[V.gene] * x[C.gene_turn]

        """
        if self.perturbation:
            for i, dv in self.perturbation.items():
                v[i] = v[i] * dv
        """

        dydt = [0] * V.NUM

        dydt[V.TGFb] = -v[1]
        dydt[V.Rec] = -v[1]
        dydt[V.TGFb_pRec] = +v[1] - v[2]
        dydt[V.S2] = -v[9] + v[11]
        dydt[V.S3] = -v[12] + v[14]
        dydt[V.S4] = (
            -3 * v[7]
            + 3 * v[8]
            - v[18]
            + v[19]
            - v[23]
            + v[24]
            - 2 * v[25]
            + 2 * v[26]
            - 2 * v[27]
            + 2 * v[28]
            - v[29]
            + v[30]
            + v[31]
        )
        dydt[V.ppS2_ppS2_ppS2] = +v[3] - v[4]
        dydt[V.ppS3_ppS3_ppS3] = +v[5] - v[6]
        dydt[V.S4_S4_S4] = +v[7] - v[8]
        dydt[V.pS2] = +v[4] + v[10] - v[11] + v[16] + v[19] + v[21] + v[26] + v[30]
        dydt[V.pS3] = +v[6] + v[13] - v[14] + v[17] + v[22] + v[24] + v[28] + v[31]
        dydt[V.ppS2] = (
            -3 * v[3]
            + 2 * v[4]
            + v[9]
            - v[10]
            - 2 * v[15]
            + v[16]
            + 2 * v[17]
            - 2 * v[18]
            + v[19]
            - v[20]
            + v[22]
            - v[25]
            - v[29]
            + v[31]
        )
        dydt[V.ppS3] = (
            -3 * v[5]
            + 2 * v[6]
            + v[12]
            - v[13]
            - v[15]
            + v[16]
            - 2 * v[20]
            + 2 * v[21]
            + v[22]
            - 2 * v[23]
            + v[24]
            - v[27]
            - v[29]
            + v[30]
        )
        dydt[V.ppS2_ppS2_S4] = +v[18] - v[19]
        dydt[V.ppS2_ppS2_ppS3] = +v[15] - v[16] - v[17]
        dydt[V.ppS2_ppS3_ppS3] = +v[20] - v[21] - v[22]
        dydt[V.ppS3_ppS3_S4] = +v[23] - v[24]
        dydt[V.ppS2_ppS3_S4] = +v[29] - v[30] - v[31]
        dydt[V.ppS3_S4_S4] = +v[27] - v[28]
        dydt[V.ppS2_S4_S4] = +v[25] - v[26]
        dydt[V.gene] = +v[32] - v[33]

        return dydt


def param_values():
    x = [0] * C.NUM

    x[C.Rec_act] = 6.00e-3
    x[C.S2tot] = 1.42e2
    x[C.S3tot] = 1.42e1
    x[C.S4tot] = 6.94e1
    x[C.S_dephos] = 2.89e-1
    x[C.S_dephosphos] = 4.95e-2
    x[C.S_phos] = 3.80e-1
    x[C.Ski_act3] = 1.41e-2
    x[C.Ski_act2] = 8.06e-4
    x[C.Ski_act1] = 0  # 1e-5
    x[C.Ski_inh3] = 0  # 1e-5
    x[C.Ski_inh2] = 4.66e-2
    x[C.Ski_inh1] = 2.68e-2
    x[C.Ski_turn] = 4.56e-3
    x[C.Skil_act3] = 4.63e-1
    x[C.Skil_act2] = 7.45e-2
    x[C.Skil_act1] = 3.07e-2
    x[C.Skil_inh3] = 0  # 1e-5
    x[C.Skil_inh2] = 7.53e-1
    x[C.Skil_inh1] = 4.10e-1
    x[C.Skil_turn] = 1.12e-2
    x[C.Dnmt3a_act3] = 7.15e-3
    x[C.Dnmt3a_act2] = 0  # 1e-5
    x[C.Dnmt3a_act1] = 0  # 1e-5
    x[C.Dnmt3a_inh3] = 0  # 1e-5
    x[C.Dnmt3a_inh2] = 1.93e-2
    x[C.Dnmt3a_inh1] = 3.62e-2
    x[C.Dnmt3a_turn] = 3.81e-3
    x[C.Sox4_act3] = 1.12e-2
    x[C.Sox4_act2] = 2.72e-4
    x[C.Sox4_act1] = 0  # 1e-5
    x[C.Sox4_inh3] = 0  # 1e-5
    x[C.Sox4_inh2] = 8.73e-2
    x[C.Sox4_inh1] = 1.03e-1
    x[C.Sox4_turn] = 8.02e-4
    x[C.Jun_act3] = 0  # 1e-5
    x[C.Jun_act2] = 1.02e0
    x[C.Jun_act1] = 5.99e0
    x[C.Jun_inh3] = 8.20e0
    x[C.Jun_inh2] = 1.46e0
    x[C.Jun_inh1] = 9.76e0
    x[C.Jun_turn] = 1.21e-1
    x[C.Smad7_act3] = 9.98e2
    x[C.Smad7_act2] = 2.12e-1
    x[C.Smad7_act1] = 2.09e1
    x[C.Smad7_inh3] = 3.67e0
    x[C.Smad7_inh2] = 0  # 1e-5
    x[C.Smad7_inh1] = 0  # 1e-5
    x[C.Smad7_turn] = 3.64e1
    x[C.Klf10_act3] = 1.00e3
    x[C.Klf10_act2] = 9.05e1
    x[C.Klf10_act1] = 0  # 1e-5
    x[C.Klf10_inh3] = 0  # 1e-5
    x[C.Klf10_inh2] = 1.79e-2
    x[C.Klf10_inh1] = 7.59e-4
    x[C.Klf10_turn] = 4.74e2
    x[C.Bmp4_act3] = 0  # 1e-5
    x[C.Bmp4_act2] = 9.04e1
    x[C.Bmp4_act1] = 1.00e3
    x[C.Bmp4_inh3] = 1.00e3
    x[C.Bmp4_inh2] = 1.94e1
    x[C.Bmp4_inh1] = 2.18e2
    x[C.Bmp4_turn] = 9.91e-1
    x[C.Cxcl15_act3] = 9.24e0
    x[C.Cxcl15_act2] = 0  # 1e-5
    x[C.Cxcl15_act1] = 0  # 1e-5
    x[C.Cxcl15_inh3] = 0  # 1e-5
    x[C.Cxcl15_inh2] = 3.20e2
    x[C.Cxcl15_inh1] = 1.00e3
    x[C.Cxcl15_turn] = 3.28e-3
    x[C.Dusp5_act3] = 0  # 1e-5
    x[C.Dusp5_act2] = 0  # 1e-5
    x[C.Dusp5_act1] = 0  # 1e-5
    x[C.Dusp5_inh3] = 4.72e-1
    x[C.Dusp5_inh2] = 1.15e-2
    x[C.Dusp5_inh1] = 1.09e-2
    x[C.Dusp5_turn] = 1.41e-2
    x[C.Tgfa_act3] = 1.28e-2
    x[C.Tgfa_act2] = 2.93e-4
    x[C.Tgfa_act1] = 0  # 1e-5
    x[C.Tgfa_inh3] = 0  # 1e-5
    x[C.Tgfa_inh2] = 1.39e0
    x[C.Tgfa_inh1] = 2.49e0
    x[C.Tgfa_turn] = 3.60e-3
    x[C.Pdk4_act3] = 0  # 1e-5
    x[C.Pdk4_act2] = 4.86e-4
    x[C.Pdk4_act1] = 0  # 1e-5
    x[C.Pdk4_inh3] = 1.29e0
    x[C.Pdk4_inh2] = 6.43e-2
    x[C.Pdk4_inh1] = 1.40e-1
    x[C.Pdk4_turn] = 1.50e-2
    x[C.k_on_223] = 0.00e0
    x[C.k_on_224] = 0.00e0
    x[C.k_on_233] = 1.47e-1
    x[C.k_on_234] = 4.03e-4
    x[C.k_on_244] = 8.63e-6
    x[C.k_on_334] = 0.00e0
    x[C.k_on_344] = 0.00e0
    x[C.k_on_222] = 0.00e0
    x[C.k_on_333] = 0.00e0
    x[C.k_on_444] = 0.00e0
    x[C.kdiss_SS] = 0.00e0
    x[C.pRec_degind] = 4.01e-2
    x[C.sd_Bmp4] = 1.41e-1
    x[C.sd_Cxcl15] = 1.35e-1
    x[C.sd_Dnmt3a] = 9.68e-2
    x[C.sd_Dusp5] = 1.06e-1
    x[C.sd_Jun] = 1.65e-1
    x[C.sd_Klf10] = 1.18e-1
    x[C.sd_Pdk4] = 7.40e-2
    x[C.sd_Ski] = 1.02e-1
    x[C.sd_Skil] = 1.75e-1
    x[C.sd_Smad7] = 7.04e-2
    x[C.sd_Sox4] = 1.26e-1
    x[C.sd_Tgfa] = 7.72e-2

    return x


def initial_values():
    y0 = [0] * V.NUM

    y0[V.TGFb] = 1
    y0[V.Rec] = 1.84
    y0[V.TGFb_pRec] = 0
    y0[V.S2] = 1.43e2
    y0[V.S3] = 1.63e1
    y0[V.S4] = 6.71e1
    y0[V.ppS2_ppS2_ppS2] = 0
    y0[V.ppS3_ppS3_ppS3] = 0
    y0[V.S4_S4_S4] = 0
    y0[V.pS2] = 0
    y0[V.pS3] = 0
    y0[V.ppS2] = 0
    y0[V.ppS3] = 0
    y0[V.ppS2_ppS2_S4] = 0
    y0[V.ppS2_ppS2_ppS3] = 0
    y0[V.ppS2_ppS3_ppS3] = 0
    y0[V.ppS3_ppS3_S4] = 0
    y0[V.ppS2_ppS3_S4] = 0
    y0[V.ppS3_S4_S4] = 0
    y0[V.ppS2_S4_S4] = 0
    y0[V.gene] = 1

    return y0