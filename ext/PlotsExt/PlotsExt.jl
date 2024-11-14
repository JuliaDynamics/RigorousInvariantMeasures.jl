module PlotsExt

using RigorousInvariantMeasures, Plots

using RecipesBase
#COV_EXCL_START
@recipe f(::Type{ApproxInducedLSV}, D::ApproxInducedLSV) = x -> plottable(D, x)
#COV_EXCL_STOP

#COV_EXCL_START
@userplot PlotC2
@recipe function f(h::PlotC2)
    if length(h.args) != 2 ||
       (typeof(h.args[1]) != C2) ||
       !(typeof(h.args[2]) <: AbstractVector)
        error("Plot C2 needs as an input a C2 Basis and a vector")
    end

    B = h.args[1]
    w = h.args[2]

    seriestype := :path
    collect(B), mid.(w)
end
#COV_EXCL_STOP

#COV_EXCL_START
using RecipesBase

"""
Plots a function in the Hat basis
"""
@recipe function f(B::Hat, w::AbstractVector)

    legend --> :bottomright

    if eltype(w) <: Interval
        w = mid.(w)
    end

    @series begin
        seriestype --> :path
        label --> L"f_{\delta}"
        ylims --> (0, NaN)
        B.p, vcat(w, w[end])
    end
end

"""
Displays error on a function in the Hat basis
"""
@recipe function f(B::Hat, error::Number, w)

    if eltype(w) <: Interval
        w = mid.(w)
    end

    if isfinite(error)
        @series begin
            seriestype --> :path
            seriesalpha --> 0.5
            fillrange --> vcat(w, w[end]) .- error
            label --> "Error area"
            B.p, vcat(w, w[end]) .+ error
        end
    end
end
#COV_EXCL_STOP

#COV_EXCL_START
using RecipesBase
using LaTeXStrings

"""
Plots a function in the Ulam basis
"""
@recipe function f(B::Ulam, w::AbstractVector)

    legend --> :bottomright

    if eltype(w) <: Interval
        w = mid.(w)
    end

    @series begin
        seriestype --> :steppost
        label --> L"f_{\delta}"
        ylims --> (0, NaN)
        B.p, vcat(w, w[end])
    end
end

"""
Displays error on a function in the Ulam basis

The w argument is unused, but kept for compatibility with other functions
for different bases
"""
@recipe function f(B::Ulam, error::Number, w = nothing)

    if isfinite(error)
        @series begin
            seriestype --> :path
            seriesalpha --> 0.5
            fillrange --> 0
            label --> "L1 Error"
            [0; sqrt(error)], [sqrt(error); sqrt(error)]
        end
    end
end
#COV_EXCL_START

using RecipesBase
#COV_EXCL_START
@recipe f(::Type{PM}, D::PM) where {PM<:PwMap} = x -> plottable(D, x)
#COV_EXCL_STOP

#COV_EXCL_START
using RecipesBase

"""
Plots a function in the Hat basis
"""
@recipe function f(B::HatNP, w::AbstractVector)

    legend --> :bottomright

    if eltype(w) <: Interval
        w = mid.(w)
    end

    @series begin
        seriestype --> :path
        label --> L"f_{\delta}"
        ylims --> (0, NaN)
        B.p, vcat(w, w[end])
    end
end

"""
Displays error on a function in the Hat basis
"""
@recipe function f(B::HatNP, error::Number, w)

    if eltype(w) <: Interval
        w = mid.(w)
    end

    if isfinite(error)
        @series begin
            seriestype --> :path
            seriesalpha --> 0.5
            fillrange --> vcat(w, w[end]) .- error
            label --> "Error area"
            B.p, vcat(w, w[end]) .+ error
        end
    end
end
#COV_EXCL_STOP

end