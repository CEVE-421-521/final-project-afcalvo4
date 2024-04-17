"""The capacity of the parking garage is 200 spaces per level"""
calculate_capacity(x::ParkingGarageState) = 200 * x.n_levels

"""
The construction will cost \$16, 000/ space for precast construction, with a 10% increase for every level above the ground level.
Note that in the adaptive case, we have to pay an additional 5% for the initial construction cost to cover the footers and columns.
"""
function calc_construction_cost(n_levels, Δn_levels, is_adaptive)
    cost_per_space = Δn_levels == 0 ? 0 : 16_000 * (1 + 0.1 * (n_levels + Δn_levels - 1))
    if is_adaptive
        cost_per_space = cost_per_space * 1.05
    end
    return cost_per_space * 200 # 200 spaces per level
end

"""Throw an error if you try to use the abstract policy type"""
function get_action(x::ParkingGarageState, policy::AbstractPolicy)
    throw("You need to implement a concrete policy type!")
end

"""
If we are following the static (fixed) policy, then the rule is simple. We add `n_levels` in the first year, and then keep the same number of levels in all subsequent years.
"""
function get_action(x::ParkingGarageState, policy::StaticPolicy)
    if x.year == 1
        return ParkingGarageAction(policy.n_levels)
    else
        return ParkingGarageAction(0)
    end
end

"""
If we are following the adaptive policy, then the rule is slightly more complicated.
We add `n_levels` in the first year. Then, every future year we compare the capacity and demand. If the demand is greater than the capacity, we add a level.
"""
function get_action(x::ParkingGarageState, policy::AdaptivePolicy)
    if x.year == 1
        return ParkingGarageAction(policy.n_levels_init)
    else
        capacity = calculate_capacity(x)
        if x.demand > capacity
            return ParkingGarageAction(1)
        else
            return ParkingGarageAction(0)
        end
    end
end

"""
Run the simulation for a single year
"""
function run_timestep(
    x::ParkingGarageState, s::ParkingGarageSOW, policy::T
) where {T<:AbstractPolicy}

    # calculate the demand for this year
    x.demand = calculate_demand(x.year, s.demand_growth_rate)

    # the very first step is to decide on the action
    a = get_action(x, policy)

    # next we have to do stuff to implement the action in our model!    
    # construction costs are higher for initial construction if we are using the adaptive policy
    is_adaptive = policy isa AdaptivePolicy && x.year == 1
    constr_costs = calc_construction_cost(x.n_levels, a.Δn_levels, is_adaptive)

    # and update the system to reflect our action!
    x.n_levels += a.Δn_levels

    # revenue -- you can only sell parking spaces that you have AND that are wanted
    capacity = calculate_capacity(x)
    revenue = 11_000 * min(capacity, x.demand)

    # lease costs are fixed
    lease_cost = 3_600_000

    # operating costs depend only on the capacity
    operating_costs = 2_000 * capacity

    # total costs
    costs = constr_costs + lease_cost + operating_costs

    profit = revenue - costs
    return profit
end

"""
Run the simulation for `s.n_years` years
"""
function simulate(s::ParkingGarageSOW, policy::T) where {T<:AbstractPolicy}

    # initialize the model
    x = ParkingGarageState()

    # calculate the profits for each year
    years = collect(1:(s.n_years))
    profits = map(years) do year
        x.year = year
        run_timestep(x, s, policy)
    end

    # apply discounting
    discount_weights = @. (1 - s.discount_rate)^(years - 1)
    npv_profits = sum(profits .* discount_weights)
    npv_profits_millions = npv_profits / 1e6
    return npv_profits_millions
end
