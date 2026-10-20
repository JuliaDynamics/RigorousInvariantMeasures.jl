@testset "infinite_sum_norms" begin
    import RigorousInvariantMeasures: infinite_sum_norms

    # Happy path: geometric tail summed from the last entry, which is < 1.
    v = vcat([1.13, 1.113, 1.047, 0.9245, 0.7349, 0.4477], fill(0.4, 23), [0.3])
    S = infinite_sum_norms(v)
    @test isfinite(S)
    @test S > 0
    # It is an upper bound for the partial sum it is built from.
    @test S >= 1.0 + sum(v[1:end-1])

    # Non-convergent data: the last norm does not certify a geometric tail, so no
    # finite bound is available.
    @test infinite_sum_norms(vcat(v[1:end-1], [1.0])) == Inf
    @test infinite_sum_norms(vcat(v[1:end-1], [1.7])) == Inf

    # Non-finite entries anywhere: Inf out, never a finite number.  The last-entry
    # case is a regression test: it used to return -0.0, which silently reported a
    # zero resolvent sum and produced spuriously tight enclosures downstream.
    vlast = copy(v); vlast[end] = Inf
    @test infinite_sum_norms(vlast) == Inf
    @test !(infinite_sum_norms(vlast) == 0.0)

    vmid = copy(v); vmid[10] = Inf
    @test infinite_sum_norms(vmid) == Inf

    vnan = copy(v); vnan[5] = NaN
    @test infinite_sum_norms(vnan) == Inf

    # Whatever the failure, the result is never a finite value on data that does
    # not support one.
    for bad in (Inf, NaN, 1.0, 2.5)
        w = copy(v); w[end] = bad
        @test !isfinite(infinite_sum_norms(w))
    end
end
