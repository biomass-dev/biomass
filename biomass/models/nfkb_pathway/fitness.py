import numpy as np
from scipy.spatial.distance import cosine

from .observable import observables, ExperimentalData, NumericalSimulation
from .set_search_param import SearchParam


def _compute_objval_rss(sim_data, exp_data):
    """Return Residual Sum of Squares"""
    return np.dot((sim_data - exp_data), (sim_data - exp_data))


def _compute_objval_cos(sim_data, exp_data):
    """Return Cosine distance"""
    return cosine(sim_data, exp_data)


def _diff_sim_and_exp(sim_matrix, exp_dict, exp_timepoint, conditions, sim_norm_max):
    sim_val = []
    exp_val = []

    for idx, condition in enumerate(conditions):
        if condition in exp_dict.keys():
            sim_val.extend(sim_matrix[list(map(int, exp_timepoint)), idx])
            exp_val.extend(exp_dict[condition])

    return np.array(sim_val) / sim_norm_max, np.array(exp_val)


def objective(indiv_gene, *args):
    """Define an objective function to be minimized"""
    if len(args) == 0:
        sp = SearchParam()
        indiv = sp.gene2val(indiv_gene)
        (x, y0) = sp.update(indiv)
    elif len(args) == 1:
        raise ValueError("not enough values to unpack (expected 2, got 1)")
    elif len(args) == 2:
        (x, y0) = args
    else:
        raise ValueError("too many values to unpack (expected 2)")

    sim = NumericalSimulation()
    exp = ExperimentalData()

    exp.set_data()

    if sim.simulate(x, y0) is None:
        error = np.zeros(len(observables))
        for i, obs_name in enumerate(observables):
            if exp.experiments[i] is not None:
                error[i] = _compute_objval_rss(
                    *_diff_sim_and_exp(
                        sim.simulations[i],
                        exp.experiments[i],
                        exp.get_timepoint(obs_name),
                        sim.conditions,
                        sim_norm_max=1
                        if not sim.normalization
                        else (
                            np.max(
                                sim.simulations[
                                    observables.index(obs_name),
                                    sim.normalization[obs_name]["timepoint"],
                                    [
                                        sim.conditions.index(c)
                                        for c in (
                                            sim.normalization[obs_name]["condition"]
                                            if sim.normalization[obs_name]["condition"]
                                            else sim.conditions
                                        )
                                    ],
                                ]
                            )
                            if sim.normalization[obs_name]["timepoint"] is not None
                            else np.max(
                                sim.simulations[
                                    observables.index(obs_name),
                                    :,
                                    [
                                        sim.conditions.index(c)
                                        for c in (
                                            sim.normalization[obs_name]["condition"]
                                            if sim.normalization[obs_name]["condition"]
                                            else sim.conditions
                                        )
                                    ],
                                ]
                            )
                        ),
                    )
                )
        return np.sum(error)  # < 1e12
    else:
        return 1e12