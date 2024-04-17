"""
A simple five-parameter model of sea-level rise (SLR).

See equation 6 of Oddo et al.

Oddo, P. C., Lee, B. S., Garner, G. G., Srikrishnan, V., Reed, P. M., Forest, C. E., & Keller, K. (2017). Deep uncertainties in sea-level rise and storm surge projections: implications for coastal flood risk management. Risk Analysis, 0(0). https://doi.org/10/ghkp82
"""

# we declare a const to tell the compiler that this value never changes, which speeds things up
# 1 mm = 0.00328084 ft
const mm_to_ft = 0.00328084

struct Oddo17SLR{T<:AbstractFloat}
    a::T
    b::T
    c::T
    tstar::T
    cstar::T
end

"""Implements Oddo equation 6 to return local sea-level in feet"""
function (s::Oddo17SLR)(t::Real)
    slr_mm =
        s.a +
        s.b * (t - 2000) +
        s.c * (t - 2000)^2 +
        s.cstar * (t > s.tstar) * (t - s.tstar)
    return slr_mm * mm_to_ft
end
