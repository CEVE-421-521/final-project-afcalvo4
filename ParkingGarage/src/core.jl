@kwdef struct ParkingGarageSOW{F<:AbstractFloat,I<:Int}
    demand_growth_rate::F = 80.0
    n_years::I = 20
    discount_rate::F = 0.12
end

"""
Demand on opening day is for 750 spaces, and rises linearly at the rate of `demand_growth_rate` spaces/ year
"""
function calculate_demand(t, demand_growth_rate::AbstractFloat)
    return 750 + demand_growth_rate * (t - 1)
end

mutable struct ParkingGarageState{T<:AbstractFloat}
    n_levels::Int
    year::Int
    demand::T
end

function ParkingGarageState()
    return ParkingGarageState(0, 1, calculate_demand(1, 80.0))
end

struct ParkingGarageAction
    Δn_levels::Int
end

abstract type AbstractPolicy end
struct StaticPolicy <: AbstractPolicy
    n_levels::Int
end
struct AdaptivePolicy <: AbstractPolicy
    n_levels_init::Int
end