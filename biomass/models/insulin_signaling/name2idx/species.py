NAMES = [
    "Ins",
    "pro_IRcom",
    "IRcom",
    "p1IRcom",
    "p2IRcom",
    "p1p2IRcom",
    "iAKT",
    "pAKT",
    "imTOR",
    "pmTOR",
    "iX",
    "pX",
    "iS6K",
    "pS6K",
    "iGSK3B",
    "pGSK3B",
    "iFoxO1",
    "pFoxO1",
    "G6Pase",
]

for idx, name in enumerate(NAMES):
    exec("{} = {:d}".format(name, idx))

NUM = len(NAMES)
