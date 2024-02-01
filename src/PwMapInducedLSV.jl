# We implement here the InducedLSV map using the new PwDynamics interface

LSV_left(x; α) = x * (1 + (2 * x)^α)
LSV_right(x) = 2 * x - 1

function PwDynamicApproxInducedLSV(α, k; T = Float64, rescaled = false)
    right = Interval{T}(1.0)

    branches = MonotonicBranch[]

    f_left(x) = LSV_left(x; α)
    f_right = LSV_right

    f = f_right
    for i = 1:k
        left = RigorousInvariantMeasures.ShootingLSV(i, 0.5, α; T = T)[1]
        push!(
            branches,
            MonotonicBranch(f, (left, right), (Interval{}(0.5), Interval{T}(1.0))),
        )
        f = f_left ∘ f
        right = left
    end

    f = x -> 0.5 / (right - 0.5) * (x - 0.5) + 0.5
    push!(
        branches,
        MonotonicBranch(f, (Interval{T}(0.5), right), (Interval{T}(0.5), Interval{T}(1.0))),
    )
    if !rescaled
        @warn "This is the induced map on [0.5, 1.0]"
    end

    return PwMap(reverse(branches), full_branch = true)
end

CoordinateChangePwMap(; T = Float64) = PwMap(
    [
        MonotonicBranch(
            x -> 2 * x - 1,
            (Interval{T}(0.5), Interval{T}(1.0)),
            (Interval{T}(0.0), Interval{T}(1.0)),
        )
    ];
    full_branch = true,
)
InvCoordinateChangePwMap(; T = Float64) = PwMap(
    [
        MonotonicBranch(
            x -> x / 2 + 0.5,
            (Interval{T}(0.0), Interval{T}(1.0)),
            (Interval{T}(0.5), Interval{T}(1.0)),
        )
    ];
    full_branch = true,
)

function RescaledApproxInducedLSV(α, k; T = Float64)
    ϕ = CoordinateChangePwMap(; T = T)
    ψ = InvCoordinateChangePwMap(; T = T)
    D = PwDynamicApproxInducedLSV(α, k; T = T, rescaled = true)
    return ϕ ∘ D ∘ ψ
end
