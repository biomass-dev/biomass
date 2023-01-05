"""
Model template for conversion of biomass models into
BioMASS.jl (https://github.com/biomass-dev/BioMASS.jl) format.

This supports BioMASS.jl v0.8.0 or more.
"""

from typing import Final

PARAMETERS: Final[
    str
] = """\
module C

const NAMES = []

for (idx,name) in enumerate(NAMES)
    eval(Meta.parse("const $name = $idx"))
end

const NUM = length(NAMES)

end  # module
"""


SPECIES: Final[
    str
] = """\
module V

const NAMES = []

for (idx, name) in enumerate(NAMES)
    eval(Meta.parse("const $name = $idx"))
end

const NUM = length(NAMES)

end  # module
"""


ODE: Final[
    str
] = """\
function diffeq!(du, u, p, t)
    v = Dict{Int64,Float64}()

    for i in 1:V.NUM
        @inbounds du[i] = 0.0
    end

end


function param_values()::Vector{Float64}
    p::Vector{Float64} = ones(C.NUM)

    return p
end


function initial_values()::Vector{Float64}
    u0::Vector{Float64} = zeros(V.NUM)

    return u0
end
"""


OBSERVABLE: Final[
    str
] = """\
const observables = []

function observables_index(observable_name::String)::Int
    if !(observable_name in observables)
        error("$observable_name is not defined in observables.")
    end
    return findfirst(isequal(observable_name),observables)
end
"""


SIMULATION: Final[
    str
] = """\
module Sim
include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./ode.jl")
include("./observable.jl")

using .C
using .V

using Sundials
using SteadyStateDiffEq

# Options for ODE solver
const ABSTOL = 1e-8
const RELTOL = 1e-8

normalization = Dict{String,Dict{}}()

const dt = 1.0
const t = collect(0.0:dt:100.0)

const conditions = []

simulations = Array{Float64,3}(
    undef, length(observables), length(conditions), length(t)
)


function solveode(
    f::Function,
    u0::Vector{Float64},
    t::Vector{Float64},
    p::Vector{Float64})::Union{ODESolution{},Nothing}
    local sol::ODESolution{}, is_successful::Bool
    prob = ODEProblem(f, u0, (t[1], t[end]), p)
    try
        prob = ODEProblem(f, u0, (t[1], t[end]), p)
        sol = solve(
            prob, CVODE_BDF(),
            abstol=ABSTOL,
            reltol=RELTOL,
            saveat=dt,
            dtmin=eps(),
            verbose=false
        )
        is_successful = ifelse(sol.t[end] == t[end], true, false)
    catch
        is_successful = false
    end
    return is_successful ? sol : nothing
end


function get_steady_state(
    f::Function,
    u0::Vector{Float64},
    p::Vector{Float64})::Vector{Float64}
    local sol::SteadyStateSolution{}
    try
        prob = ODEProblem(f, u0, (0.0, Inf), p)
        prob = SteadyStateProblem(prob)
        sol = solve(
            prob,
            DynamicSS(
                CVODE_BDF();
                abstol=ABSTOL,
                reltol=RELTOL
            ),
            dt=dt,
            dtmin=eps(),
            verbose=false
        )
        #is_successful = ifelse(sol.retcode === :Success, true, false)
        return sol.u
    catch
        return []
    end
    #return is_successful ? sol.u : []
end


function simulate!(p::Vector{Float64}, u0::Vector{Float64})::Union{Bool,Nothing}
    # unperturbed steady state

    # add ligand
    for (i, condition) in enumerate(conditions)

        sol = solveode(diffeq!, u0, t, p)
        if sol === nothing
            return false
        else
            @inbounds @simd for j in eachindex(t)
                # line_num + 4
            end
        end
    end
end
end # module
"""


EXPERIMENTAL_DATA: Final[
    str
] = """\
module Exp
include("./observable.jl")

experiments = Array{Dict{String,Array{Float64,1}},1}(undef, length(observables))
error_bars = Array{Dict{String,Array{Float64,1}},1}(undef, length(observables))


function get_timepoint(obs_name::String)::Vector{Float64}
    return []
end
end # module
"""


SEARCH_PARAM: Final[
    str
] = """\
# Specify model parameters and/or initial values to optimize
function get_search_index()::Tuple{Array{Int64,1},Array{Int64,1}}
    # parameters
    search_idx_params::Vector{Int} = []

    # initial values
    search_idx_initials::Vector{Int} = []

    return search_idx_params, search_idx_initials
end


function get_search_region()::Matrix{Float64}
    p::Vector{Float64} = param_values()
    u0::Vector{Float64} = initial_values()

    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = get_search_index()
    search_param::Vector{Float64} = initialize_search_param(search_idx, p, u0)

    search_rgn::Matrix{Float64} = zeros(2, length(p) + length(u0))

    # Default: 0.1 ~ 10x
    for (i, j) in enumerate(search_idx[1])
        search_rgn[1,j] = search_param[i] * 0.1  # lower bound
        search_rgn[2,j] = search_param[i] * 10.0  # upper bound
    end

    # Default: 0.5 ~ 2x
    for (i, j) in enumerate(search_idx[2])
        search_rgn[1,j + length(p)] = search_param[i + length(search_idx[1])] * 0.5  # lower bound
        search_rgn[2,j + length(p)] = search_param[i + length(search_idx[1])] * 2.0  # upper bound
    end

    # search_rgn[:, C.param_name] = [lower_bound, upper_bound]
    # search_rgn[:, V.var_name+length(p)] = [lower_bound, upper_bound]

    search_rgn = convert_scale(search_rgn, search_idx)

    return search_rgn
end


function update_param(indiv::Vector{Float64})::Tuple{Array{Float64,1},Array{Float64,1}}
    p::Vector{Float64} = param_values()
    u0::Vector{Float64} = initial_values()

    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = get_search_index()

    for (i, j) in enumerate(search_idx[1])
        @inbounds p[j] = indiv[i]
    end
    for (i, j) in enumerate(search_idx[2])
        @inbounds u0[j] = indiv[i + length(search_idx[1])]
    end

    # parameter constraints

    return p, u0
end


function decode_gene2val(indiv_gene)::Vector{Float64}
    search_rgn::Matrix{Float64} = get_search_region()
    indiv::Vector{Float64} = zeros(length(indiv_gene))

    for i in eachindex(indiv_gene)
        indiv[i] = 10^(
            indiv_gene[i] * (
                search_rgn[2,i] - search_rgn[1,i]
            ) + search_rgn[1,i]
        )
    end

    return round.(indiv, sigdigits=7)
end


function encode_val2gene(indiv::Vector{Float64})
    search_rgn::Matrix{Float64} = get_search_region()
    indiv_gene::Vector{Float64} = zeros(length(indiv))

    for i in eachindex(indiv)
        indiv_gene[i] = (
            log10(indiv[i]) - search_rgn[1,i]
        ) / (
            search_rgn[2,i] - search_rgn[1,i]
        )
    end

    return indiv_gene
end


function encode_bestIndivVal2randGene(
        gene_idx::Int64,
        best_indiv::Vector{Float64},
        p0_bounds::Vector{Float64})::Float64
    search_rgn::Matrix{Float64} = get_search_region()
    rand_gene::Float64 = (
        log10(
            best_indiv[gene_idx] * 10^(
                rand() * log10(p0_bounds[2] / p0_bounds[1]) + log10(p0_bounds[1])
            )
        ) - search_rgn[1,gene_idx]
    ) / (
        search_rgn[2,gene_idx] - search_rgn[1,gene_idx]
    )
    return rand_gene
end


function initialize_search_param(
        search_idx::Tuple{Array{Int64,1},Array{Int64,1}},
        p::Vector{Float64},
        u0::Vector{Float64})::Vector{Float64}
    duplicate::Vector{String} = []
    if length(search_idx[1]) != length(unique(search_idx[1]))
        for idx in findall(
            [count(x -> x == i, search_idx[1]) for i in unique(search_idx[1])] .!= 1
        )
            push!(duplicate, C.NAMES[search_idx[1][idx]])
        end
        error(
            "Duplicate parameters (C.): $duplicate"
        )
    elseif length(search_idx[2]) != length(unique(search_idx[2]))
        for idx in findall(
            [count(x -> x == i, search_idx[2]) for i in unique(search_idx[2])] .!= 1
        )
            push!(duplicate, V.NAMES[search_idx[2][idx]])
        end
        error(
            "Duplicate initial conditions (V.): $duplicate"
        )
    end
    search_param = zeros(
        length(search_idx[1]) + length(search_idx[2])
    )
    for (i, j) in enumerate(search_idx[1])
        @inbounds search_param[i] = p[j]
    end
    for (i, j) in enumerate(search_idx[2])
        @inbounds search_param[i + length(search_idx[1])] = u0[j]
    end

    if any(x -> x == 0.0, search_param)
        msg::String = "search_param must not contain zero."
        for idx in search_idx[1]
            if p[idx] == 0.0
                error(
                    @sprintf(
                        "`C.%s` in search_idx_params: ", C.NAMES[idx]
                    ) * msg
                )
            end
        end
        for idx in search_idx[2]
            if u0[idx] == 0.0
                error(
                    @sprintf(
                        "`V.%s` in search_idx_initials: ", V.NAMES[idx]
                    ) * msg
                )
            end
        end
    end

    return search_param
end


function convert_scale(
        search_rgn::Matrix{Float64},
        search_idx::Tuple{Array{Int64,1},Array{Int64,1}})::Matrix{Float64}
    for i = 1:size(search_rgn, 2)
        if minimum(search_rgn[:,i]) < 0.0
            msg = "search_rgn[lower_bound,upper_bound] must be positive.\n"
            if i <= C.NUM
                error(@sprintf("`C.%s` ", C.NAMES[i]) * msg)
            else
                error(@sprintf("`V.%s` ", V.NAMES[i - C.NUM]) * msg)
            end
        elseif minimum(search_rgn[:,i]) == 0.0 && maximum(search_rgn[:,i]) != 0.0
            msg = "lower_bound must be larger than 0.\n"
            if i <= C.NUM
                error(@sprintf("`C.%s` ", C.NAMES[i]) * msg)
            else
                error(@sprintf("`V.%s` ", V.NAMES[i - C.NUM]) * msg)
            end
        elseif search_rgn[2,i] - search_rgn[1,i] < 0.0
            msg = "lower_bound must be smaller than upper_bound.\n"
            if i <= C.NUM
                error(@sprintf("`C.%s` ", C.NAMES[i]) * msg)
            else
                error(@sprintf("`V.%s` ", V.NAMES[i - C.NUM]) * msg)
            end
        end
    end

    nonzero_idx::Vector{Int} = []
    for i = 1:size(search_rgn, 2)
        if search_rgn[:,i] != [0.0,0.0]
            push!(nonzero_idx, i)
        end
    end
    difference::Vector{Int} = collect(
        symdiff(
            Set(nonzero_idx),
            Set(append!(search_idx[1], C.NUM .+ search_idx[2]))
        )
    )
    if length(difference) > 0
        for idx in difference
            if idx <= C.NUM
                println(@sprintf("`C.%s`", C.NAMES[Int(idx)]))
            else
                println(@sprintf("`V.%s`", V.NAMES[Int(idx) - C.NUM]))
            end
        end
        error(
            "Set these search_params in both search_idx and search_rgn."
        )
    end

    search_rgn = search_rgn[:,nonzero_idx]

    return log10.(search_rgn)
end
"""

PROBLEM: Final[
    str
] = """\
# Residual Sum of Squares
function compute_objval_rss(
    sim_data::Vector{Float64},
    exp_data::Vector{Float64})::Float64
    error::Float64 = 0.0
    for i in eachindex(exp_data)
        @inbounds error += (sim_data[i] - exp_data[i])^2
    end
    return error
end


# Cosine similarity
function compute_objval_cos(
    sim_data::Vector{Float64},
    exp_data::Vector{Float64})::Float64
    error::Float64 = 1.0 - dot(sim_data, exp_data) / (norm(sim_data) * norm(exp_data))
    return error
end


function conditions_index(condition_name::String)::Int
    if !(condition_name in Sim.conditions)
        error("$condition_name is not defined in Sim.conditions")
    end
    return findfirst(isequal(condition_name), Sim.conditions)
end


function diff_sim_and_exp(
    sim_matrix::Matrix{Float64},
    exp_dict::Dict{String,Array{Float64,1}},
    exp_timepoint::Vector{Float64},
    conditions::Vector{String};
    sim_norm_max::Float64)::Tuple{Vector{Float64},Vector{Float64}}
    sim_result::Vector{Float64} = []
    exp_result::Vector{Float64} = []

    for (idx, condition) in enumerate(conditions)
        if condition in keys(exp_dict)
            append!(sim_result, sim_matrix[idx, Int.(exp_timepoint .+ 1)])
            append!(exp_result, exp_dict[condition])
        end
    end

    return (sim_result ./ sim_norm_max, exp_result)
end


# Define an objective function to be minimized.
function objective(indiv_gene)::Float64
    indiv::Vector{Float64} = decode_gene2val(indiv_gene)

    (p, u0) = update_param(indiv)

    if Sim.simulate!(p, u0) isa Nothing
        error::Vector{Float64} = zeros(length(observables))
        for (i, obs_name) in enumerate(observables)
            if isassigned(Exp.experiments, i)
                if length(Sim.normalization) > 0
                    norm_max::Float64 = (
                        Sim.normalization[obs_name]["timepoint"] !== nothing ? maximum(
                            Sim.simulations[
                                i,
                                [conditions_index(c) for c in Sim.normalization[obs_name]["condition"]],
                                Sim.normalization[obs_name]["timepoint"]
                            ]
                        ) : maximum(
                            Sim.simulations[
                                i,
                                [conditions_index(c) for c in Sim.normalization[obs_name]["condition"]],
                                :,
                            ]
                        )
                    )
                end
                error[i] = compute_objval_rss(
                    diff_sim_and_exp(
                        Sim.simulations[i, :, :],
                        Exp.experiments[i],
                        Exp.get_timepoint(obs_name),
                        Sim.conditions,
                        sim_norm_max=ifelse(
                            length(Sim.normalization) == 0, 1.0, norm_max
                        )
                    )...
                )
            end
        end
        return sum(error) # < 1e12
    else
        return 1e12
    end
end
"""
