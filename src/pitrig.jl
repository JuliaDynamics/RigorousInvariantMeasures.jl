"""
Modifies the implementation of sinpi() and cospi() from IntervalArithmetic
to improve it so that sinpi(1) == 0
"""

using ValidatedNumerics

using ValidatedNumerics.IntervalArithmetic: atomic, SVector, @round

function find_quadrantspi(x::T) where {T}
    temp = IntervalArithmetic.atomic(Interval{T}, x) * 2

    return IntervalArithmetic.SVector(floor(temp.lo), floor(temp.hi))
end

sinpi(a...) = Base.Math.sinpi(a...)

function sinpi(a::Interval{T}) where T
    isempty(a) && return a

    whole_range = Interval{T}(-1, 1)

    diam(a) > 2 && return whole_range

    # The following is equiavlent to doing temp = a / half_pi  and
    # taking floor(a.lo), floor(a.hi)
    lo_quadrant = minimum(find_quadrantspi(a.lo))
    hi_quadrant = maximum(find_quadrantspi(a.hi))

    if hi_quadrant - lo_quadrant > 4  # close to limits
        return whole_range
    end

    lo_quadrant = mod(lo_quadrant, 4)
    hi_quadrant = mod(hi_quadrant, 4)

    # Different cases depending on the two quadrants:
    if lo_quadrant == hi_quadrant
        a.hi - a.lo > 1 && return whole_range  # in same quadrant but separated by almost 2pi
        lo = @round(sinpi(a.lo), sinpi(a.lo)) # Interval(sin(a.lo, RoundDown), sin(a.lo, RoundUp))
        hi = @round(sinpi(a.hi), sinpi(a.hi)) # Interval(sin(a.hi, RoundDown), sin(a.hi, RoundUp))
        return hull(lo, hi)

    elseif lo_quadrant==3 && hi_quadrant==0
        return @round(sinpi(a.lo), sinpi(a.hi)) # Interval(sin(a.lo, RoundDown), sin(a.hi, RoundUp))

    elseif lo_quadrant==1 && hi_quadrant==2
        return @round(sinpi(a.hi), sinpi(a.lo)) # Interval(sin(a.hi, RoundDown), sin(a.lo, RoundUp))

    elseif ( lo_quadrant == 0 || lo_quadrant==3 ) && ( hi_quadrant==1 || hi_quadrant==2 )
        return @round(min(sinpi(a.lo), sinpi(a.hi)), 1)
        # Interval(min(sin(a.lo, RoundDown), sin(a.hi, RoundDown)), one(T))

    elseif ( lo_quadrant == 1 || lo_quadrant==2 ) && ( hi_quadrant==3 || hi_quadrant==0 )
        return @round(-1, max(sinpi(a.lo), sinpi(a.hi)))
        # Interval(-one(T), max(sin(a.lo, RoundUp), sin(a.hi, RoundUp)))

    else #if( lo_quadrant == 0 && hi_quadrant==3 ) || ( lo_quadrant == 2 && hi_quadrant==1 )
        return whole_range
    end
end

cospi(a...) = Base.Math.cospi(a...)

function cospi(a::Interval{T}) where T
    isempty(a) && return a

    whole_range = Interval(-one(T), one(T))

    diam(a) > 2 && return whole_range

    lo_quadrant = minimum(find_quadrantspi(a.lo))
    hi_quadrant = maximum(find_quadrantspi(a.hi))

    if hi_quadrant - lo_quadrant > 4  # close to limits
        return Interval(-one(T), one(T))
    end

    lo_quadrant = mod(lo_quadrant, 4)
    hi_quadrant = mod(hi_quadrant, 4)

    # Different cases depending on the two quadrants:
    if lo_quadrant == hi_quadrant # Interval limits in the same quadrant
        a.hi - a.lo > 1 && return whole_range
        lo = @round(cospi(a.lo), cospi(a.lo))
        hi = @round(cospi(a.hi), cospi(a.hi))
        return hull(lo, hi)

    elseif lo_quadrant == 2 && hi_quadrant==3
        return @round(cospi(a.lo), cospi(a.hi))

    elseif lo_quadrant == 0 && hi_quadrant==1
        return @round(cospi(a.hi), cospi(a.lo))

    elseif ( lo_quadrant == 2 || lo_quadrant==3 ) && ( hi_quadrant==0 || hi_quadrant==1 )
        return @round(min(cospi(a.lo), cospi(a.hi)), 1)

    elseif ( lo_quadrant == 0 || lo_quadrant==1 ) && ( hi_quadrant==2 || hi_quadrant==3 )
        return @round(-1, max(cospi(a.lo), cospi(a.hi)))

    else#if ( lo_quadrant == 3 && hi_quadrant==2 ) || ( lo_quadrant == 1 && hi_quadrant==0 )
        return whole_range
    end
end
