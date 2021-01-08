from .name2idx import C, V


class DifferentialEquation(object):
    def __init__(self, perturbation):
        self.perturbation = perturbation

    def diffeq(self, t, y, *x):
        """Table S6: List of differential equations of the reduced model"""
        dydt = [0] * V.NUM

        dydt[V.TNFR] = 1 * x[C.uptake] * x[C.TNF] - 1 * x[C.deact_TNFR] * y[V.TNFR]
        dydt[V.Ikk] = -1 * (x[C.act_Ikk_by_TNF] * y[V.TNFR] * y[V.Ikk]) + 1 * (x[C.trigger_iIkk] * y[V.iIkk])
        dydt[V.pIkk] = 1 * (x[C.act_Ikk_by_TNF] * y[V.TNFR] * y[V.Ikk]) - 1 * (x[C.act_pIkk] * y[V.pIkk])
        dydt[V.ppIkk] = 1 * (x[C.act_pIkk] * y[V.pIkk]) - 1 * (x[C.deact_ppIkk] * y[V.ppIkk])
        dydt[V.iIkk] = 1 * (x[C.deact_ppIkk] * y[V.ppIkk]) - 1 * (x[C.trigger_iIkk] * y[V.iIkk])
        dydt[V.NfkIkb] = (
            -1 * ((x[C.act_Ikb_by_Ikk]) * y[V.pIkk] * y[V.NfkIkb])
            - 1 * ((x[C.act_Nfk_by_Ikk]) * y[V.pIkk] * y[V.NfkIkb])
            + 1 * (x[C.form_complex] * y[V.Nfk] * y[V.Ikb])
            + 1 * (x[C.ext_nNfkIkb] * y[V.nNfkIkb]) * (x[C.Vnuc] / 1)
        )
        dydt[V.NfkpIkb] = (
            1 * ((x[C.act_Ikb_by_Ikk]) * y[V.pIkk] * y[V.NfkIkb])
            - 1 * ((x[C.act_Nfk_by_Ikk_complex]) * y[V.pIkk] * y[V.NfkpIkb])
            - 1 * (x[C.split_NfkpIkb] * y[V.NfkpIkb])
        )
        dydt[V.pNfkIkb] = (
            -1 * ((x[C.act_Ikb_complex]) * y[V.pNfkIkb])
            - 1 * ((x[C.act_Ikb_by_Ikk]) * y[V.pIkk] * y[V.pNfkIkb])
            + 1 * ((x[C.act_Nfk_by_Ikk]) * y[V.pIkk] * y[V.NfkIkb])
        )
        dydt[V.pNfkpIkb] = (
            1 * ((x[C.act_Ikb_complex]) * y[V.pNfkIkb])
            + 1 * ((x[C.act_Ikb_by_Ikk]) * y[V.pIkk] * y[V.pNfkIkb])
            + 1 * ((x[C.act_Nfk_by_Ikk_complex]) * y[V.pIkk] * y[V.NfkpIkb])
            - 1 * (x[C.split_NfkIkb] * y[V.pNfkpIkb])
        )
        dydt[V.pNfk] = 1 * (x[C.split_NfkIkb] * y[V.pNfkpIkb]) - 1 * ((x[C.int_Nfk] * x[C.eta_int_pNfk]) * y[V.pNfk])
        dydt[V.Nfk] = (
            1 * (x[C.split_NfkpIkb] * y[V.NfkpIkb])
            - 1 * (x[C.form_complex] * y[V.Nfk] * y[V.Ikb])
            - 1 * ((x[C.int_Nfk]) * y[V.Nfk])
        )
        dydt[V.pIkb] = (
            1 * (x[C.split_NfkpIkb] * y[V.NfkpIkb])
            + 1 * (x[C.split_NfkIkb] * y[V.pNfkpIkb])
            - 1 * (x[C.degrad_Ikb] * y[V.pIkb])
        )
        dydt[V.Ikb] = (
            -1 * (x[C.form_complex] * y[V.Nfk] * y[V.Ikb])
            + 1 * (x[C.prod_Ikb] * y[V.mIkb])
            - 1 * (x[C.int_Ikb] * y[V.Ikb])
        )
        dydt[V.mIkb] = 1 * (x[C.prod_mIkb_by_nNfk] * y[V.nNfk]) - 1 * (x[C.degrad_mIkb] * y[V.mIkb])
        dydt[V.nIkb] = 1 * (x[C.int_Ikb] * y[V.Ikb]) * (1 / x[C.Vnuc]) - 1 * (
            x[C.form_complex_nuc] * y[V.nNfk] * y[V.nIkb]
        )
        dydt[V.pnNfk] = 1 * ((x[C.int_Nfk] * x[C.eta_int_pNfk]) * y[V.pNfk]) * (1 / x[C.Vnuc]) - 1 * (
            x[C.deact_pnNfk] * y[V.pnNfk]
        )
        dydt[V.nNfk] = (
            1 * ((x[C.int_Nfk]) * y[V.Nfk]) * (1 / x[C.Vnuc])
            + 1 * (x[C.deact_pnNfk] * y[V.pnNfk])
            - 1 * (x[C.form_complex_nuc] * y[V.nNfk] * y[V.nIkb])
        )
        dydt[V.nNfkIkb] = 1 * (x[C.form_complex_nuc] * y[V.nNfk] * y[V.nIkb]) - 1 * (x[C.ext_nNfkIkb] * y[V.nNfkIkb])
        dydt[V.RnaA20_1] = 1 * (x[C.build_RnaA20] * y[V.nNfk]) - 1 * (x[C.shuttle_RnaA20] * y[V.RnaA20_1])
        dydt[V.RnaA20] = 1 * (x[C.shuttle_RnaA20] * y[V.RnaA20_1]) - 1 * (x[C.degrad_RnaA20] * y[V.RnaA20])
        dydt[V.A20] = 1 * (x[C.build_A20] * y[V.RnaA20]) - 1 * (x[C.degrad_A20] * y[V.A20])

        return dydt


def param_values():
    x = [0] * C.NUM

    x[C.uptake] = 1.0000
    x[C.TNF] = 1.0000
    x[C.trigger_iIkk] = 0.0041
    x[C.deact_TNFR] = 0.0010
    x[C.deact_ppIkk] = 0.1660
    x[C.deact_pnNfk] = 1000.0000
    x[C.act_Ikk_by_TNF] = 0.0714
    x[C.act_pIkk] = 0.0648
    x[C.act_Ikb_by_Ikk] = 0.3980
    x[C.act_Nfk_by_Ikk] = 0.6438
    x[C.act_Nfk_by_Ikk_complex] = 0.2816
    x[C.act_Ikb_complex] = 1.3897
    x[C.form_complex] = 2.8390
    x[C.form_complex_nuc] = 1000.0000
    x[C.ext_nNfkIkb] = 1000.0000
    x[C.Vnuc] = 1.0000
    x[C.split_NfkpIkb] = 0.0811
    x[C.split_NfkIkb] = 1.0000
    x[C.int_Nfk] = 0.0100
    x[C.int_Ikb] = 0.1226
    x[C.eta_int_pNfk] = 17.9585
    x[C.degrad_Ikb] = 0.6308
    x[C.degrad_mIkb] = 0.0313
    x[C.degrad_RnaA20] = 0.0089
    x[C.degrad_A20] = 0.0116
    x[C.prod_Ikb] = 1.0000
    x[C.prod_mIkb_by_nNfk] = 0.0047
    x[C.build_RnaA20] = 1.0000
    x[C.build_A20] = 0.0006
    x[C.shuttle_RnaA20] = 0.0311

    return x


def initial_values():
    y0 = [0] * V.NUM

    y0[V.Ikk] = 1.0
    y0[V.NfkIkb] = 1.0

    return y0