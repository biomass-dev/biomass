NAMES = [
    'MP',
    'MC',
    'MB',
    'PC',
    'CC',
    'PCP',
    'CCP',
    'PCC',
    'PCN',
    'PCCP',
    'PCNP',
    'BC',
    'BCP',
    'BN',
    'BNP',
    'IN',
]

for idx, name in enumerate(NAMES):
    exec(
        '{} = {:d}'.format(
            name, idx
        )
    )

NUM = len(NAMES)