NAMES = [
    "uptake",
    "TNF",
    "trigger_iIkk",
    "deact_TNFR",
    "deact_ppIkk",
    "deact_pnNfk",
    "act_Ikk_by_TNF",
    "act_pIkk",
    "act_Ikb_by_Ikk",
    "act_Nfk_by_Ikk",
    "act_Nfk_by_Ikk_complex",
    "act_Ikb_complex",
    "form_complex",
    "form_complex_nuc",
    "ext_nNfkIkb",
    "Vnuc",
    "split_NfkpIkb",
    "split_NfkIkb",
    "int_Nfk",
    "int_Ikb",
    "eta_int_pNfk",
    "degrad_Ikb",
    "degrad_mIkb",
    "degrad_RnaA20",
    "degrad_A20",
    "prod_Ikb",
    "prod_mIkb_by_nNfk",
    "build_RnaA20",
    "build_A20",
    "shuttle_RnaA20",
]

for idx, name in enumerate(NAMES):
    exec("{} = {:d}".format(name, idx))

NUM = len(NAMES)