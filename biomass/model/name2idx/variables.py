var_names = [
    'ppMEKc',
    'CREBn',
    'pCREBn',
    'ERKc',
    'ERKn',
    'pERKc',
    'pERKn',
    'ppERKc',
    'ppERKn',
    'Elk1n',
    'pElk1n',
    'cFOSc',
    'cFOSn',
    'pcFOSc',
    'pcFOSn',
    'DUSPc',
    'DUSPn',
    'pDUSPc',
    'pDUSPn',
    'DUSPn_ERKn',
    'DUSPn_pERKn',
    'DUSPn_ppERKn',
    'pDUSPn_ERKn',
    'pDUSPn_pERKn',
    'pDUSPn_ppERKn',
    'RSKc',
    'pRSKc',
    'pRSKn',
    'PrecfosmRNAn',
    'PreduspmRNAn',
    'cfosmRNAc',
    'duspmRNAc',
    'Fc',
    'Fn',
    'FmRNAc',
    'PreFmRNAn',
]

for idx, name in enumerate(var_names):
    exec(
        '{} = {:d}'.format(
            name, idx
        )
    )

len_f_vars = len(var_names)