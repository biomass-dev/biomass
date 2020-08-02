NAMES = [
    'MKKK',
    'MKKK_P',
    'MKK',
    'MKK_P',
    'MKK_PP',
    'MAPK',
    'MAPK_P',
    'MAPK_PP',
]

for idx, name in enumerate(NAMES):
    exec(
        '{} = {:d}'.format(
            name, idx
        )
    )

NUM = len(NAMES)