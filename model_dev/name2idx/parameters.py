"""
This is where you define the names of the parameters.

NAMES : list of strings

"""
NAMES = [

]

for idx, name in enumerate(NAMES):
    exec(
        '{} = {:d}'.format(
            name, idx
        )
    )

NUM = len(NAMES)
