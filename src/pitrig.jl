"""
Modifies the implementation of sinpi() and cospi() from IntervalArithmetic
to improve it so that sinpi(1) == 0
"""



using IntervalArithmetic: atomic, SVector, @round

function find_quadrantspi(x::T) where {T}
    temp = IntervalArithmetic.atomic(Interval{T}, x) * 2

    return IntervalArithmetic.SVector(floor(inf(temp)), floor(sup(temp)))
end

sinpi(a...) = Base.Math.sinpi(a...)

function sinpi(a::Interval{T}) where {T}
    isempty_interval(a) && return a

    whole_range = interval(T, -1, 1)

    diam(a) > 2 && return whole_range

    # The following is equiavlent to doing temp = a / half_pi  and
    # taking floor(inf(a)), floor(sup(a))
    lo_quadrant = minimum(find_quadrantspi(inf(a)))
    hi_quadrant = maximum(find_quadrantspi(sup(a)))

    if hi_quadrant - lo_quadrant > 4  # close to limits
        return whole_range
    end

    lo_quadrant = mod(lo_quadrant, 4)
    hi_quadrant = mod(hi_quadrant, 4)

    # Different cases depending on the two quadrants:
    if lo_quadrant == hi_quadrant
        sup(a) - inf(a) > 1 && return whole_range  # in same quadrant but separated by almost 2pi
        lo = @round(T, sinpi(inf(a)), sinpi(inf(a))) # interval(sin(inf(a), RoundDown), sin(inf(a), RoundUp))
        hi = @round(T, sinpi(sup(a)), sinpi(sup(a))) # interval(sin(sup(a), RoundDown), sin(sup(a), RoundUp))
        return hull(lo, hi)

    elseif lo_quadrant == 3 && hi_quadrant == 0
        return @round(T, sinpi(inf(a)), sinpi(sup(a))) # interval(sin(inf(a), RoundDown), sin(sup(a), RoundUp))

    elseif lo_quadrant == 1 && hi_quadrant == 2
        return @round(T, sinpi(sup(a)), sinpi(inf(a))) # interval(sin(sup(a), RoundDown), sin(inf(a), RoundUp))

    elseif (lo_quadrant == 0 || lo_quadrant == 3) && (hi_quadrant == 1 || hi_quadrant == 2)
        return @round(T, min(sinpi(inf(a)), sinpi(sup(a))), 1)
        # interval(min(sin(inf(a), RoundDown), sin(sup(a), RoundDown)), one(T))

    elseif (lo_quadrant == 1 || lo_quadrant == 2) && (hi_quadrant == 3 || hi_quadrant == 0)
        return @round(T, -1, max(sinpi(inf(a)), sinpi(sup(a))))
        # interval(-one(T), max(sin(inf(a), RoundUp), sin(sup(a), RoundUp)))

    else #if( lo_quadrant == 0 && hi_quadrant==3 ) || ( lo_quadrant == 2 && hi_quadrant==1 )
        return whole_range
    end
end

cospi(a...) = Base.Math.cospi(a...)

function cospi(a::Interval{T}) where {T}
    isempty_interval(a) && return a

    whole_range = interval(-one(T), one(T))

    diam(a) > 2 && return whole_range

    lo_quadrant = minimum(find_quadrantspi(inf(a)))
    hi_quadrant = maximum(find_quadrantspi(sup(a)))

    if hi_quadrant - lo_quadrant > 4  # close to limits
        return interval(-one(T), one(T))
    end

    lo_quadrant = mod(lo_quadrant, 4)
    hi_quadrant = mod(hi_quadrant, 4)

    # Different cases depending on the two quadrants:
    if lo_quadrant == hi_quadrant # Interval limits in the same quadrant
        sup(a) - inf(a) > 1 && return whole_range
        lo = @round(T, cospi(inf(a)), cospi(inf(a)))
        hi = @round(T, cospi(sup(a)), cospi(sup(a)))
        return hull(lo, hi)

    elseif lo_quadrant == 2 && hi_quadrant == 3
        return @round(T, cospi(inf(a)), cospi(sup(a)))

    elseif lo_quadrant == 0 && hi_quadrant == 1
        return @round(T, cospi(sup(a)), cospi(inf(a)))

    elseif (lo_quadrant == 2 || lo_quadrant == 3) && (hi_quadrant == 0 || hi_quadrant == 1)
        return @round(T, min(cospi(inf(a)), cospi(sup(a))), 1)

    elseif (lo_quadrant == 0 || lo_quadrant == 1) && (hi_quadrant == 2 || hi_quadrant == 3)
        return @round(T, -1, max(cospi(inf(a)), cospi(sup(a))))

    else#if ( lo_quadrant == 3 && hi_quadrant==2 ) || ( lo_quadrant == 1 && hi_quadrant==0 )
        return whole_range
    end
end

# These definitions were missing, too

sinpi(x::Taylor1) = sin(x * π)
cospi(x::Taylor1) = cos(x * π)
