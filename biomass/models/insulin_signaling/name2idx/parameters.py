NAMES = [
    "k1_synthesis",
    "k1_InsIRcom",
    "k2_InsIRcom",
    "k1_p1IRcomDeg",
    "k1_p1IRcomPhos",
    "k1_p1p2IRcomdePhos",
    "k1_IRcomPhos",
    "k1_p2IRcomdePhos",
    "k1_p2IRcomDeg",
    "k1_Insp2IRcom",
    "k2_Insp2IRcom",
    "k1_p1p2IRcomDeg",
    "k1_AKTPhos",
    "k1_pAKTdePhos",
    "k1_mTORPhos",
    "k1_pmTORdePhos",
    "k1_S6KPhos",
    "k1_pS6KdePhos",
    "k1_XPhos",
    "k1_pXdePhos",
    "k1_GSK3BPhos",
    "k1_pGSK3BdePhos",
    "k1_FoxO1Phos",
    "k1_pFoxO1dePhos",
    "k1_G6PaseSynthesis",
    "k1_G6PaseDeg",
]

for idx, name in enumerate(NAMES):
    exec("{} = {:d}".format(name, idx))

NUM = len(NAMES)
