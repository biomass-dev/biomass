NAMES = [
    'V1',
    'n',
    'KI',
    'K1',
    'V2',
    'K2',
    'k3',
    'K3',
    'k4',
    'K4',
    'V5',
    'K5',
    'V6',
    'K6',
    'k7',
    'K7',
    'k8',
    'K8',
    'V9',
    'K9',
    'V10',
    'K10',
]

for idx, name in enumerate(NAMES):
    exec(
        '{} = {:d}'.format(
            name, idx
        )
    )

NUM = len(NAMES)
