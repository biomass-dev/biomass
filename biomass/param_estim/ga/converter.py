import numpy as np


def decode_gene2variable(individual_gene, search_rgn):
    indiv_var = 10**(
        individual_gene * (
            search_rgn[1, :] - search_rgn[0, :]
        ) + search_rgn[0, :]
    )

    return indiv_var


def encode_bestindiv2randgene(best_indiv, search_rgn, p0_bounds):
    rand_gene = (
        np.log10(
            best_indiv * 10**(
                np.random.rand(len(best_indiv))
                * np.log10(p0_bounds[1]/p0_bounds[0])
                + np.log10(p0_bounds[0])
            )
        ) - search_rgn[0, :]
    ) / (search_rgn[1, :] - search_rgn[0, :])

    return rand_gene
