using Base: @kwdef
using DataFrames
using Interpolations
using Unitful

struct DepthDamageFunction{I<:Interpolations.AbstractExtrapolation}
    # the interpolation object
    itp::I
end

function DepthDamageFunction(depths_ft::Vector{T}, damages::Vector{T}) where {T<:Real}
    itp = Interpolations.LinearInterpolation(
        depths_ft, damages; extrapolation_bc=Interpolations.Flat()
    )
    return DepthDamageFunction(itp)
end

"""A house has a value in USD and an area"""
@kwdef struct House{T<:Real}

    # essential features
    area_ft2::T
    value_usd::T
    height_above_gauge_ft::T

    # depth-damage
    ddf::DepthDamageFunction

    # some metadata fields
    occupancy::AbstractString
    dmg_fn_id::AbstractString
    source::AbstractString
    description::AbstractString
    comment::AbstractString
end

"""
Constructor function to parse a row from a DataFrame, eg from `haz_fl_dept.csv` and convert it to a `House` object
"""
function House(
    row::DataFrames.DataFrameRow;
    area::Unitful.Area,
    height_above_gauge::Unitful.Length,
    value_usd::T,
) where {T<:Real}

    # get metadata fields
    occupancy = string(row.Occupancy)
    dmg_fn_id = string(row.DmgFnId)
    source = string(row.Source)
    description = string(row.Description)
    comment = string(row.Comment)

    # now get the depths and damages
    depths_ft = Float64[]
    damages = Float64[]

    # Iterate over each column in the row
    for (col_name, value) in pairs(row)
        # Check if the column name starts with 'ft'
        if startswith(string(col_name), "ft")
            if value != "NA"
                # Convert column name to depth
                depth_str = string(col_name)[3:end] # Extract the part after 'ft'
                is_negative = endswith(depth_str, "m")
                depth_str = is_negative ? replace(depth_str, "m" => "") : depth_str
                depth_str = replace(depth_str, "_" => ".") # Replace underscore with decimal
                depth_val = parse(Float64, depth_str)
                depth_val = is_negative ? -depth_val : depth_val
                push!(depths_ft, depth_val)

                # Add value to damages
                push!(damages, parse(Float64, string(value)))
            end
        end
    end

    ddf = DepthDamageFunction(depths_ft, damages)

    return House(;
        area_ft2=ustrip(u"ft^2", area),
        height_above_gauge_ft=ustrip(u"ft", height_above_gauge),
        value_usd=value_usd,
        ddf=ddf,
        occupancy=occupancy,
        dmg_fn_id=dmg_fn_id,
        source=source,
        description=description,
        comment=comment,
    )
end

"""Call a depth-damage function on a house"""
function (ddf::DepthDamageFunction)(depth_ft::AbstractFloat)
    return ddf.itp(depth_ft)
end

function (ddf::DepthDamageFunction)(depth::Unitful.Length)
    depth_ft = ustrip(u"ft", depth)
    return ddf(depth_ft)
end

"""
A struct to define the cost of elevating a house
    
See Zarekarizi et al. (2020) or Doss-Gollin and Keller (2023) for details

Zarekarizi, M., Srikrishnan, V., & Keller, K. (2020). Neglecting uncertainties biases house-elevation decisions to manage riverine flood risks. Nature Communications, 11(1), 5361. https://doi.org/10.1038/s41467-020-19188-9
Doss-Gollin, J., & Keller, K. (2023). A subjective Bayesian framework for synthesizing deep uncertainties in climate risk management. Earth’s Future, 11(1). https://doi.org/10.1029/2022EF003044
"""
struct ElevationCostCalculator
    # at its heart it's just an interpolator!
    itp::Interpolations.Extrapolation
end

# Constructor to initialize the struct with default values and set up the interpolation
function ElevationCostCalculator()
    elevation_thresholds = [0.0, 5.0, 8.5, 12.0, 14.0]
    elevation_rates = [80.36, 82.5, 86.25, 103.75, 113.75]
    itp = LinearInterpolation(elevation_thresholds, elevation_rates)
    return ElevationCostCalculator(itp)
end

# Define the method for calculating the elevation cost
function (calculator::ElevationCostCalculator)(house::House, Δh_ft::T) where {T<:Real}

    # cannot lower the house
    Δh_ft < 0.0 && throw(DomainError(Δh, "Cannot lower the house"))

    # cannot elevate >14ft
    Δh_ft > 14.0 && throw(DomainError(Δh, "Cannot elevate >14ft"))

    # no cost if no elevation
    Δh_ft ≈ 0.0 && return 0.0

    # otherwise, we're going to calculate the cost
    base_cost = (10000 + 300 + 470 + 4300 + 2175 + 3500) # in USD
    rate = calculator.itp(Δh_ft)
    cost = base_cost + house.area_ft2 * rate
    return cost
end

# Define the method for calculating the elevation cost
function (calculator::ElevationCostCalculator)(
    house::House, Δh::T
) where {T<:Unitful.Length}

    # Convert Δh to feet
    Δh_ft = ustrip(u"ft", Δh)

    # Call the method for calculating the elevation cost
    return calculator(house, Δh_ft)
end

elevation_cost = ElevationCostCalculator()