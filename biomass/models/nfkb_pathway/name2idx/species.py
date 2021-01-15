NAMES = [
    "TNFR",
    "Ikk",
    "pIkk",
    "ppIkk",
    "iIkk",
    "NfkIkb",
    "NfkpIkb",
    "pNfkIkb",
    "pNfkpIkb",
    "pNfk",
    "Nfk",
    "pIkb",
    "Ikb",
    "mIkb",
    "nIkb",
    "pnNfk",
    "nNfk",
    "nNfkIkb",
    "RnaA20_1",
    "RnaA20",
    "A20",
]

for idx, name in enumerate(NAMES):
    exec("{} = {:d}".format(name, idx))

NUM = len(NAMES)