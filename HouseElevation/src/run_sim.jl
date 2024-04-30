using Distributions

"""Helper function for trapezoidal rule"""
function trapz(x, y)
    return sum((x[2:end] - x[1:(end - 1)]) .* (y[2:end] + y[1:(end - 1)])) * 0.5
end

"""
Run the model for a given action and SOW

Expected Annual Damages are computed using the trapezoidal rule
"""
function run_sim(a::Action, sow::SOW, p::ModelParams)

    # first, we calculate the cost of elevating the house
    construction_cost = elevation_cost(p.house, a.Δh_ft)

    # we don't need to recalculate the steps of the trapezoidal integral for each year
    storm_surges_ft = range(
        quantile(sow.surge_dist, 0.0005); stop=quantile(sow.surge_dist, 0.9995), length=130
    )

    eads = map(p.years) do year

        # get the sea level for this year
        slr_ft = sow.slr(year)

        # Compute EAD using trapezoidal rule
        pdf_values = pdf.(sow.surge_dist, storm_surges_ft) # probability of each
        depth_ft_gauge = storm_surges_ft .+ slr_ft # flood at gauge
        depth_ft_house = depth_ft_gauge .- (p.house.height_above_gauge_ft + a.Δh_ft) # flood @ house
        damages_frac = p.house.ddf.(depth_ft_house) ./ 100 # damage
        weighted_damages = damages_frac .* pdf_values # weighted damage
        # Trapezoidal integration of weighted damages
        ead = trapz(storm_surges_ft, weighted_damages) * p.house.value_usd
    end

    years_idx = p.years .- minimum(p.years)
    discount_fracs = (1 - sow.discount_rate) .^ years_idx
    ead_npv = sum(eads .* discount_fracs)
    return -(ead_npv + construction_cost)
end

# The proactive function triggers an elevation action only if the freboard distance Δ is surpass
function get_ProactiveAction(x::BuildingState, policy::ProactivePolicy, a::Action)
    if (x.elevation - x.MSL) < policy.delta_min
        return BuildingAction(a.Δh_ft)
    else
        return BuildingAction(0)
    end
    
end

# The proactive function triggers an elevation action only if the maximum level threshold Lmax is surpass
function get_ReactiveAction(x::BuildingState, policy::ReactivePolicy, a::Action)
    if x.max_level > policy.critical_level
        return BuildingAction(a.Δh_ft)
    else
        return BuildingAction(0)
    end
    
end

# One time step in the sequential decision siluation
function run_timestep(
    x::BuildingState, sow::SOW, a::Action, p::ModelParams, policy::T
) where {T<:AbstractPolicy}

    # Initial level of the building
    if x.year == 1
        x.elevation = p.house.height_above_gauge_ft
    end

    # Mean sea level from the SLR model
    x.MSL = sow.slr(p.years[x.year])

    # Set the elevation action according to a policy
    if policy isa ProactivePolicy
        Δ = get_ProactiveAction(x, policy, a)
    elseif policy isa ReactivePolicy
        Δ = get_ReactiveAction(x, policy, a)
    end

    # Estimates the construction cost of the elevation action
    construction_cost = elevation_cost(p.house, Δ.Δelevation)

    # Change the elevation state of the building
    x.elevation += Δ.Δelevation

    # Calculate the maximum yearly level or year flood event
    x.max_level = rand(sow.surge_dist) + x.MSL - x.elevation - Δ.Δelevation
    # Estimates the building damages and losses
    damages_year = p.house.ddf(x.max_level)/100
    x.Damage = damages_year
    EAD = x.Damage * p.house.value_usd
    # The year cost cashflow is the sum of the construction losses costs
    cost = EAD + construction_cost
    
    return cost
end

function simulate(sow::SOW, a::Action, p::ModelParams, policy::T) where {T<:AbstractPolicy}

    # initialize the model
    x = BuildingState()
    
    # Calculate the year cashflow
    years = collect(1:size(p.years)[1])
    costs = map(years) do year
        x.year = year
        run_timestep(x, sow, a, p, policy)
    end

    # Calculate the net present value of the simullation
    discount_weights = @. (1 - sow.discount_rate)^(years - 1)
    npv = sum(costs .* discount_weights)
    npv_millions = npv

    return npv_millions
end
