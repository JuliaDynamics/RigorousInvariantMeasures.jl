using RigorousInvariantMeasures, IntervalArithmetic, TaylorSeries, TaylorModels, Test
const IA = IntervalArithmetic

@testset "Jets and Taylor models of the transfer operator" begin
    # encloses(J, t): every coefficient of the BigFloat truth t lies in that of J
    tobig(x::Interval) = interval(BigFloat, IA.inf(x), IA.sup(x))
    encloses(J, t, k) = all(issubset_interval(t[m], tobig(J[m])) for m = 0:k)
    k = 6

    @testset "doubling map, f = x^3: Lf(y) = (y^3 + (y+1)^3)/16" begin
        D = mod1_dynamic(x -> 2x)
        truth3(y) = let t = Taylor1([interval(BigFloat, y), interval(BigFloat, 1)], k)
            (t^3 + (t + 1)^3) / 16
        end
        J = transfer_jet(D, x -> x^3, 0.3, k)
        @test encloses(J, truth3(0.3), k)
        @test all(IA.diam(J[m]) < 1e-14 for m = 0:k)
        # a circle map reads the point modulo 1
        C = CircleMap(D, :ω)
        @test encloses(transfer_jet(C, x -> x^3, 1.3, k), truth3(0.3), k)
    end

    # T(x) = 2x + x(1-x)/2, inverse branches τ_j(y) = (5 - sqrt(25 - 8(y+j)))/2,
    # |τ_j'(y)| = 2/sqrt(25 - 8(y+j)), j = 0, 1
    D = mod1_dynamic(x -> 2x + x * (1 - x) / 2; full_branch = true)
    g(x) = exp(x) + x^2
    truth(y) = let t = Taylor1([interval(BigFloat, y), interval(BigFloat, 1)], k)
        sum(g((5 - sqrt(25 - 8(t + j))) / 2) * (2 / sqrt(25 - 8(t + j))) for j = 0:1)
    end

    @testset "jets at a point and on an interval, against the explicit branches" begin
        J = transfer_jet(D, g, 0.3, k)
        @test encloses(J, truth(0.3), k)
        @test all(IA.diam(J[m]) < 1e-13 for m = 0:k)
        JI = transfer_jet(D, g, interval(0.3, 0.32), k)
        @test all(encloses(JI, truth(y), k) for y in (0.3, 0.31, 0.32))
        JB = transfer_jet(D, g, interval(BigFloat, 0.3), k)
        # the truth at 512 bits, finer than the 256 bits of JB
        tB = setprecision(BigFloat, 512) do
            truth(0.3)
        end
        @test all(issubset_interval(tB[m], JB[m]) for m = 0:k)
        @test all(IA.diam(JB[m]) < 1e-60 for m = 0:k)
        @test_throws ArgumentError transfer_jet(D, g, 0.3, -1)
    end

    @testset "Taylor model of Lf on an interval" begin
        I = interval(0.3, 0.32)
        M = transfer_taylor_model(D, g, I, k)
        @test IA.diam(remainder(M)) < 1e-10
        for y in (0.3, 0.305, 0.31, 0.3175, 0.32)
            @test issubset_interval(truth(y)[0], tobig(M(interval(y) - interval(IA.mid(I)))))
        end
    end

    @testset "a branch that does not cover the point" begin
        # x -> 1.5x mod 1: branch 1 on [0, 2/3] onto [0, 1], branch 2 on [2/3, 1] onto
        # [0, 1/2]; Lf for f = 1 is 2/1.5 on [0, 1/2) and 1/1.5 on (1/2, 1]
        E = mod1_dynamic(x -> 1.5x)
        @test in_interval(1 / 1.5, transfer_jet(E, x -> one(x), 0.7, 2)[0])
        @test in_interval(2 / 1.5, transfer_jet(E, x -> one(x), 0.3, 2)[0])
        @test in_interval(1 / 1.5,
                          transfer_taylor_model(E, x -> one(x), interval(0.6, 0.8), 3)(interval(0.05)))
        @test_throws ArgumentError transfer_jet(E, x -> one(x), interval(0.45, 0.55), 2)
    end

    # every cell of every piece with cells encloses the truth L at up to three interior points
    function cells_enclose(Ps, L)
        ok = true
        for P in Ps, q in (P.cells === nothing ? [] : P.cells)
            y0 = expansion_point(q.model)
            a, b = IA.inf(q.piece), IA.sup(q.piece)
            for y in (a + (b - a) / 5, (a + b) / 2, b - (b - a) / 5)
                a < y < b || continue            # interior points only (tiny cells round)
                ok &= issubset_interval(L(interval(BigFloat, y)), tobig(q.model(interval(y) - y0)))
            end
        end
        return ok
    end
    ends(Ps) = [IA.sup(P.piece) for P in Ps]

    # x -> 1.5x mod 1: (Lg)(y) = (g(y/1.5) + [y < 1/2] g((y+1)/1.5)) / 1.5
    E = mod1_dynamic(x -> 1.5x)
    th = interval(BigFloat, 1.5)
    Lb(h) = y -> (h(y / th) + (y < 0.5 ? h((y + 1) / th) : interval(BigFloat, 0))) / th
    # T(x) = 2x + x(1-x)/2, explicit branches
    τ(y, j) = (5 - sqrt(25 - 8(y + j))) / 2
    dτ(y, j) = 2 / sqrt(25 - 8(y + j))
    Lt(h) = y -> sum(h(τ(y, j)) * dτ(y, j) for j = 0:1)

    @testset "pieces cut at the branch image endpoints" begin
        Ps = transfer_taylor_models(E, exp, interval(0.3, 0.7), 5)
        @test ends(Ps) == [0.5, 0.7] && all(length(P.cells) == 1 for P in Ps)
        @test cells_enclose(Ps, Lb(exp))
        @test all(IA.diam(remainder(P.cells[1].model)) < 1e-8 for P in Ps)   # measured 1.5e-9, 3.6e-10
        # x -> √2 x mod 1: the second image ends at an enclosure of √2 - 1, not thin
        s2 = sqrt(interval(2))
        Q = transfer_taylor_models(mod1_dynamic(x -> s2 * x), exp, interval(0.3, 0.6), 4)
        @test length(Q) == 3 && Q[2].cells === nothing
        @test in_interval(sqrt(big(2)) - 1, interval(BigFloat, IA.inf(Q[2].piece), IA.sup(Q[2].piece)))
        r2 = sqrt(interval(BigFloat, 2))
        Lg(two) = y -> (exp(y / r2) + (two ? exp((y + 1) / r2) : interval(BigFloat, 0))) / r2
        @test cells_enclose(Q[1:1], Lg(true)) && cells_enclose(Q[3:3], Lg(false))
    end

    @testset "pieces between cuts partitioned into cells" begin
        R1 = transfer_taylor_models(D, g, interval(0.0, 1.0), 5; parts = 8)
        @test length(R1) == 1 && length(R1[1].cells) == 8
        @test all(IA.sup(R1[1].cells[i].piece) == IA.inf(R1[1].cells[i+1].piece) for i = 1:7)
        @test cells_enclose(R1, Lt(g))
        @test maximum(IA.diam(remainder(q.model)) for q in R1[1].cells) < 1e-6   # measured 5.7e-8
        S1 = transfer_taylor_models(E, exp, interval(0.0, 1.0), 5; parts = 3)
        @test ends(S1) == [0.5, 1.0] && all(length(P.cells) == 3 for P in S1)
    end

    @testset "iterating on a piecewise f" begin
        k = 5
        # 1.5x mod 1: L(exp) jumps at 1/2, so L²(exp) also jumps at T(1/2) = 3/4
        P1 = transfer_taylor_models(E, exp, interval(0.0, 1.0), k)
        P2 = transfer_taylor_models(E, P1, interval(0.0, 1.0), k)
        @test ends(P2) == [0.5, 0.75, 1.0] && all(P.cells !== nothing for P in P2)
        @test cells_enclose(P2, Lb(Lb(exp)))
        # with cells: the cell boundaries of L(exp) give cell boundaries of L²(exp), not cuts
        C1 = transfer_taylor_models(E, exp, interval(0.0, 1.0), k; parts = 3)
        C2 = transfer_taylor_models(E, C1, interval(0.0, 1.0), k; parts = 3)
        @test ends(C2) == [0.5, 0.75, 1.0] && all(P.cells !== nothing for P in C2)
        @test cells_enclose(C2, Lb(Lb(exp)))
        # the full-branch map: no cut inside (0, 1) at any step
        Q1 = transfer_taylor_models(D, g, interval(0.0, 1.0), k; parts = 8)
        Q2 = transfer_taylor_models(D, Q1, interval(0.0, 1.0), k; parts = 8)
        Q3 = transfer_taylor_models(D, Q2, interval(0.0, 1.0), k; parts = 8)
        @test length(Q2) == 1 && length(Q3) == 1
        @test cells_enclose(Q2, Lt(Lt(g))) && cells_enclose(Q3, Lt(Lt(Lt(g))))
        @test_throws ArgumentError transfer_taylor_models(E, P1, interval(0.0, 1.0), k + 1)
    end
end
