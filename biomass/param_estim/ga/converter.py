import numpy as np

def decode_gene2variable(individual_gene,search_region):

    return 10.**(individual_gene*(search_region[1,:] - search_region[0,:]) + search_region[0,:])


def encode_bestindiv2randgene(best_indiv,search_region,p0_bounds):
    
    return (
            np.log10(
                best_indiv*10**(
                    np.random.rand(len(best_indiv))*\
                        np.log10(p0_bounds[1]/p0_bounds[0])+\
                            np.log10(p0_bounds[0])
                )
            ) - search_region[0,:]
            )/(search_region[1,:] - search_region[0,:])