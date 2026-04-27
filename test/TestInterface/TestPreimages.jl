using Test
using RigorousInvariantMeasures: preimages, is_increasing
using IntervalArithmetic

@testset "preimages" begin

    # Testing skip_beginning / skip_end

    a1 = LinRange(0, 1, 12)
    a2 = [
        0,
        interval(0, 0.1),
        interval(0, 0.2),
        interval(0.1, 0.3),
        interval(0.2, 0.4),
        interval(0.3, 0.5),
        interval(0.3, 0.6),
        interval(0.6, 0.8),
        interval(0.8, 1),
        interval(1, 1),
    ] #tricky weakly sorted vector

    x1 = interval(0.25, 0.25)
    x2 = interval(0.1, 0.8)
    x3 = 2
    x4 = -1
    x5 = 0
    x6 = 0.3

    # test for the for wide intervals and 0 derivative 

    @test isequal_interval(
        RigorousInvariantMeasures.range_estimate_monotone(x -> x, interval(0, 1)),
        interval(0, 1),
    )
    @test isequal_interval(
        RigorousInvariantMeasures.range_estimate_monotone(x -> exp(x), interval(0, 1)),
        interval(1, sup(exp(interval(1)))),
    )

    begin
        preim_branch = RigorousInvariantMeasures.preimage
        br = RigorousInvariantMeasures.MonotonicBranch(
            x -> x^2,
            (@interval(0), @interval(1)),
        )

        ϵ = 1.0 / 1024
        max_iter = 100


        # if X is the search interval we test when issubset_interval(f(X), y)
        root = preim_branch(interval(-5, 5), br, @interval(0.0, 1.0); ϵ, max_iter = 100)
        @test issubset_interval(root, (@interval(0.0, 1.0) + @interval(-ϵ, ϵ)))

        # if X is the search interval we test when issubset_interval(y, f(X))
        root = preim_branch(interval(0.0, 0.04), br, @interval(0.0, 1.0); ϵ, max_iter = 100)
        @test issubset_interval(root, (@interval(0, 0.2) + @interval(-ϵ, ϵ)))

        # if X is the search interval we test when y ∩ f(X) ≂̸ ∅, the expected result is
        # f^{-1}(y) ∩ X
        root = preim_branch(interval(0.0, 0.04), br, @interval(0.1, 1.0); ϵ, max_iter = 100)
        @test issubset_interval(root, (@interval(0.1, 0.2) + @interval(-ϵ, ϵ)))

        root = preim_branch(interval(0.0, 0.04), br, @interval(0.2, 1.0); ϵ, max_iter = 100)
        @test issubset_interval(root, (@interval(0.2) + @interval(-ϵ, ϵ)))

        # if X is the search interval we test when y ∩ f(X) = ∅ 
        root = preim_branch(interval(0.0, 0.04), br, @interval(0.3, 1.0); ϵ, max_iter = 10)
        @test isempty_interval(root)

        # we test when the domain of f contains the zero of the derivative

        br = RigorousInvariantMeasures.MonotonicBranch(
            x -> (x - @interval(0.1))^2,
            (@interval(0.1), @interval(1)),
        )

        # if X is the search interval we test when issubset_interval(y, f(X))
        root = preim_branch(interval(0.0, 0.04), br, @interval(0.1, 1.0); ϵ, max_iter = 100)
        @test issubset_interval(root, (@interval(0.1, 0.3) + @interval(-ϵ, ϵ)))

        # if X is the search interval we test when y ∩ f(X) ≂̸ ∅, the expected result is
        # f^{-1}(y) ∩ X
        root = preim_branch(interval(0.0, 0.04), br, @interval(0.2, 1.0); ϵ, max_iter = 100)
        @test issubset_interval(root, (@interval(0.2, 0.3) + @interval(-ϵ, ϵ)))

        root = preim_branch(interval(0.0, 0.04), br, @interval(0.0); ϵ, max_iter = 100)
        @test issubset_interval(root, (@interval(0.0) + @interval(-ϵ, ϵ)))

        # if X is the search interval we test when y ∩ f(X) = ∅ 
        root = preim_branch(interval(0.0, 0.04), br, @interval(0.4, 1.0); ϵ, max_iter = 10)
        @test isempty_interval(root)

        # we test the exit rule for Krawczyk, i.e., if in_interval(0, f)′(x_mid)
        # return x
        # remark that to bypass the unique increasing test we need to call 
        # the full MonotonicBranch constructor. In general this function would not be allowed
        # by the constructor. MonotonicBranch are monotone in interval(X[1].hi, X[2].lo)
        br = RigorousInvariantMeasures.MonotonicBranch(
            x -> x^2,
            (@interval(-1), @interval(1)),
            (interval(1), interval(1)),
            true,
        )
        root = preim_branch(interval(0.0), br, @interval(0.0); ϵ, max_iter = 10)
        @test isequal_interval(root, interval(0.0))


    end

    for a in (a1, a2)
        for x in (x1, x2, x3, x4, x5, x6)
            i = RigorousInvariantMeasures.first_overlapping(a, x)
            @test i == 0 || sup(interval(a[i])) <= inf(interval(x))
            @test i == length(a) || !(sup(interval(a[i+1])) <= inf(interval(x)))

            j = RigorousInvariantMeasures.last_overlapping(a, x)
            @test j == 0 || !(inf(interval(a[j])) >= sup(interval(x)))
            @test j == length(a) || inf(interval(a[j+1])) >= sup(interval(x))
        end
    end

    # checking for equality with looser tolerance
    function approxintervals(a, b)
        return inf(a) ≈ inf(b) && sup(a) ≈ sup(b)
    end

    for y in (a1, a2)
        for f in (x -> x / 2, x -> 1 + x / 2, x -> 1 - x / 2)
            b = RigorousInvariantMeasures.MonotonicBranch(f, (@interval(0), @interval(1)))
            ylabel = 1:length(y)
            x, xlabel =
                RigorousInvariantMeasures.preimages(y, b, ylabel; ϵ = 1e-14, max_iter = 100)
            @test isequal_interval(x[1], b.X[1])
            @test !isequal_interval(x[end], b.X[2]) # to make sure the last entry isn't there
            @test length(x) == length(xlabel)
            z1 = filter(
                x -> !isempty_interval(x),
                intersect_interval.(map(interval, y), hull(b.Y[1], b.Y[2] + 1e-15)),
            ) # the 1e-15 is there to simulate "disjoint intervals"
            z2 = f.(x)
            if !is_increasing(b)
                z2 = reverse(z2)
            end

            @test all(approxintervals.(z1, z2))
        end
    end

    # Test nonmonotone 
    f = x -> (x - @interval(0.3))^2
    b = RigorousInvariantMeasures.MonotonicBranch(f, (@interval(0.3), @interval(1)))
    x = RigorousInvariantMeasures.preimage(
        0.04,
        b,
        interval(0.2, 0.7),
        ϵ = 1e-14,
        max_iter = 100,
    )
    @test in_interval(0.5, x)

    b = RigorousInvariantMeasures.MonotonicBranch(f, (@interval(0), @interval(0.3)))
    x = RigorousInvariantMeasures.preimage(
        0.04,
        b,
        interval(0, 0.4),
        ϵ = 1e-14,
        max_iter = 100,
    )
    @test in_interval(0.1, x)

    #@test_logs (:debug,"Not contracting, fallback to bisection") RigorousInvariantMeasures.preimage(0.04, b, interval(0, 0.4), ϵ = 1e-14, max_iter = 100)

    D = mod1_dynamic(x -> 2x)
    DD = ∘(D, D, D, D)
    p = [0, 0.2, 0.4, 0.6, 0.8]

    x, xlabel = preimages(p, DD; ϵ = 1e-13, max_iter = 100)
    @test all(in_interval.(0:1/80:79/80, x))
    @test xlabel == repeat([1, 2, 3, 4, 5], 16)

    # test composing two different functions, to check that composition is handled in the correct order

    f = x -> 2x
    g = x -> (x + x^2) / 2

    D1 = mod1_dynamic(f)
    D2 = mod1_dynamic(g)

    @test RigorousInvariantMeasures.is_increasing(D1)

    y = 0:0.2:1
    x, xlabel = preimages(y, D1 ∘ D2; ϵ = 1e-13, max_iter = 100)

    @test all(in_interval.(0:0.2:1.8, f.(g.(x))))
    @test xlabel == repeat(1:5, 2)

    D1 = PwMap(
        [x -> 2x, x -> 6x - 3, x -> 3x - 2],
        [0, 0.5, @interval(2 / 3), 1],
        [0 1; 0 1; 0 1],
    )
    D2 = PwMap([x -> 2x, x -> 4x - 2, x -> 4x - 3], [0, 0.5, 0.75, 1], [0 1; 0 1; 0 1])

    @test RigorousInvariantMeasures.is_increasing(D1)

    z = 0:0.3:1
    y, ylabel, y′ =
        RigorousInvariantMeasures.preimages_and_derivatives(z, D1; ϵ = 0.0, max_iter = 100)

    @test all(in_interval.(
        [
            0,
            0.15,
            0.3,
            0.45,
            0.5,
            0.5 + 0.05,
            0.5 + 0.1,
            0.5 + 0.15,
            2 / 3,
            2 / 3 + 0.1,
            2 / 3 + 0.2,
            2 / 3 + 0.3,
        ],
        y,
    ))
    @test ylabel == repeat(1:4, 3)
    @test all(in_interval.([2, 2, 2, 2, 6, 6, 6, 6, 3, 3, 3, 3], y′))

    x, xlabel, x′ = RigorousInvariantMeasures.preimages_and_derivatives(
        z,
        D1 ∘ D2;
        ϵ = 0.0,
        max_iter = 100,
    )

    @test all(in_interval.(vcat(mid.(y) / 2, 0.5 .+ mid.(y) / 4, 0.75 .+ mid.(y) / 4), x))
    @test x′ == kron([4, 12, 6, 8, 24, 12, 8, 24, 12], [1, 1, 1, 1])

    f = x -> x^2
    g = x -> 1 - x^2
    Dinc = PwMap([f], [interval(0, 0), interval(1, 1)])
    Ddec = PwMap([g], [interval(0, 0), interval(1, 1)])

    @test is_increasing(Dinc)
    @test !is_increasing(Ddec)
    @test is_increasing(Ddec ∘ Ddec)
    @test !is_increasing(Dinc ∘ Ddec)
    @test !is_increasing(Ddec ∘ Ddec ∘ Ddec)

    z = [interval(0, 0), interval(0.5, 0.5), interval(1, 1)]
    x, xlabel, x′ = RigorousInvariantMeasures.preimages_and_derivatives(
        z,
        Dinc ∘ Ddec;
        ϵ = 0.0,
        max_iter = 100,
    )

    @test all(in_interval.([0, sqrt(1 - 1 / sqrt(2))], x))
    @test all(xlabel .== [2, 1])
    @test all(in_interval.([0, -2 * sqrt(2 - sqrt(2))], x′))

    x, xlabel, x′ = RigorousInvariantMeasures.preimages_and_derivatives(
        z,
        Ddec ∘ Ddec;
        ϵ = 0.0,
        max_iter = 100,
    )
    @test all(in_interval.([0, sqrt(1 - 1 / sqrt(2))], x))
    @test all(xlabel .== [1, 2])
    @test all(in_interval.([0, 2 * sqrt(2 - sqrt(2))], x′))

    D = PwMap(
        [
            x -> 17 * x / 5,
            x -> (34 * ((17 * x - 5) / 17) / 25 + 3) * ((17 * x - 5) / 17),
            x -> (34 * ((17 * x - 10) / 17) / 25 + 3) * ((17 * x - 10) / 17),
            x -> 17 * ((17 * x - 15) / 17) / 5,
        ],
        [interval(0), interval(5) / 17, interval(10) / 17, interval(15) / 17, interval(1)],
        [
            interval(0) interval(1)
            interval(0) interval(1)
            interval(0) interval(1)
            interval(0) @interval(0.4)
        ],
    )

    # we just check that this doesn't throw, for now
    RigorousInvariantMeasures.preimages(0:0.25:1, D; ϵ = 0.0, max_iter = 100)

    # test has_infinite_derivative_at_endpoints
    using RigorousInvariantMeasures: has_infinite_derivative_at_endpoints

    D = PwMap(
        [
            x -> 17 * x / 5,
            x -> (34 * ((17 * x - 5) / 17) / 25 + 3) * ((17 * x - 5) / 17),
            x -> (34 * ((17 * x - 10) / 17) / 25 + 3) * ((17 * x - 10) / 17),
            x -> 17 * ((17 * x - 15) / 17) / 5,
        ],
        [interval(0), interval(5) / 17, interval(10) / 17, interval(15) / 17, interval(1)],
        [
            interval(0) interval(1)
            interval(0) interval(1)
            interval(0) interval(1)
            interval(0) @interval(0.4)
        ],
    )

    for i = 1:length(D.branches)
        @test has_infinite_derivative_at_endpoints(D.branches[i]) == (false, false)
    end
    @test !has_infinite_derivative_at_endpoints(D)

    D0 = Lorenz()
    D = D0 ∘ D0 ∘ D0

    @test has_infinite_derivative_at_endpoints(D.E.branches[1]) == (false, true)
    for i = 2:7
        @test has_infinite_derivative_at_endpoints(D.E.branches[i]) == (true, true)
    end
    @test has_infinite_derivative_at_endpoints(D.E.branches[end]) == (true, false)

    @test has_infinite_derivative_at_endpoints(D.E)

    D = BZ()
    # @test has_infinite_derivative_at_endpoints(D)

end #testset
