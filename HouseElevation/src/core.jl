using Base: @kwdef
using Distributions
using StatsBase: mean

"""ModelParams contains all the variables that are constant across simulations"""
@kwdef struct ModelParams
    house::House
    years::Vector{Int}
end

"""A SOW contains all the variables that may vary from one simulation to the next"""
struct SOW{T<:Real}
    slr::Oddo17SLR # the parameters of sea-level rise
    surge_dist::Distributions.UnivariateDistribution # the distribution of storm surge
    discount_rate::T # the discount rate, as a percentage (e.g., 2% is 0.02)
end

"""
In this model, we only hvae one decision variable: how high to elevate the house.
"""
struct Action{T<:Real}
    Δh_ft::T
end
function Action(Δh::T) where {T<:Unitful.Length}
    Δh_ft = ustrip(u"ft", Δh)
    return Action(Δh_ft)
end

# BuildingState x is a state variable that saves each time step information
mutable struct BuildingState{T<:AbstractFloat}
    elevation::T
    year::Int
    MSL::T
    Damage::T
    max_level::T
end

# The function is inicialized in zeros except for the year in 1
function BuildingState()
    return BuildingState(0.0, 1, 0.0, 0.0, 0.0)
end

# The building action is an elevation increment Δ
struct BuildingAction
    Δelevation::AbstractFloat
end

# Two types of abstract policies are defined
abstract type AbstractPolicy end
# The proactive policy uses a Free-board distance to trigger a decision sequence
struct ProactivePolicy <: AbstractPolicy
    delta_min::AbstractFloat
end
# The reactive policy uses a year maximum water level to trigger a decision sequence
struct ReactivePolicy <: AbstractPolicy
    critical_level::AbstractFloat
end