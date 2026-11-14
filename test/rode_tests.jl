using LagrangianDescriptors
using StochasticDiffEq
using QuadGK, Test
using SciMLBase: EnsembleSerial, remake
using LagrangianDescriptors: augmentprob
using LinearAlgebra: norm
using Statistics: cor
using DiffEqNoiseProcess: WienerProcess
using SciMLBase: EnsembleProblem

# The integral of the values `vals` at the times `t` by the trapezoidal rule, which matches how
# `RandomHeun` integrates the Lagrangian descriptors, so the comparisons can use a tight tolerance
trapezoid(vals, t) = sum((vals[begin + 1:end] .+ vals[begin:end - 1]) ./ 2 .* diff(t))

@testset "Augmented RODE" begin
    @testset "Linear OOP" begin
        f = function (u, p, t, W)
            p * u + W
        end
        p = 0.1
        t0 = 0.0
        tf = 5.0
        tspan = (t0, tf)
        u0 = 0.5
        prob = RODEProblem(f, u0, tspan, p)

        M = function (du, u, p, t, W)
            sum(abs2, du)
        end

        augprob = @test_nowarn augmentprob(prob, M)
        augsol = @test_nowarn solve(augprob, RandomHeun(), dt=1/100, save_noise=true)
        @test length(augsol.u[begin]) == 4
        @test augsol.u[begin].fwd == augsol.u[begin].bwd == u0
        @test augsol.t ≈ augsol.W.t
        @test augsol.u[end].lfwd ≈ trapezoid(M.(f.(getindex.(augsol.u, :fwd), p, augsol.t, first.(augsol.W.u)), getindex.(augsol.u, :fwd), p, augsol.t, first.(augsol.W.u)), augsol.t) rtol = 1e-3
        # the backward direction has its own noise, the last entry of the augmented noise
        @test augsol.u[end].lbwd ≈ trapezoid(M.(f.(getindex.(augsol.u, :bwd), p, augsol.t, last.(augsol.W.u)), getindex.(augsol.u, :bwd), p, augsol.t, last.(augsol.W.u)), augsol.t) rtol = 1e-3
        @test augsol.u[begin].lfwd ≈ 0.0
        @test augsol.u[begin].lbwd ≈ 0.0

        augprob = @test_nowarn augmentprob(prob, M; direction = :forward)
        augsol = @test_nowarn solve(augprob, RandomHeun(), dt=1/100, save_noise=true)
        @test length(augsol.u[begin]) == 2
        @test augsol.u[begin].fwd == u0
        @test augsol.u[begin].lfwd ≈ 0.0
        @test augsol.u[end].lfwd ≈ trapezoid(M.(f.(getindex.(augsol.u, :fwd), p, augsol.t, first.(augsol.W.u)), getindex.(augsol.u, :fwd), p, augsol.t, first.(augsol.W.u)), augsol.t) rtol = 1e-3

        augprob = @test_nowarn augmentprob(prob, M; direction = :backward)
        augsol = @test_nowarn solve(augprob, RandomHeun(), dt=1/100, save_noise=true)
        @test length(augsol.u[begin]) == 2
        @test augsol.u[begin].bwd == u0
        @test augsol.u[begin].lbwd ≈ 0.0
        # the backward direction has its own noise, the last entry of the augmented noise
        @test augsol.u[end].lbwd ≈ trapezoid(M.(f.(getindex.(augsol.u, :bwd), p, augsol.t, last.(augsol.W.u)), getindex.(augsol.u, :bwd), p, augsol.t, last.(augsol.W.u)), augsol.t) rtol = 1e-3

        @test_throws ArgumentError augmentprob(prob, M; direction = :blah)
    end

    @testset "Cubic IIP" begin
        f = function (u)
            u .- u .^ 3
        end
        f! = function (du, u, p, t, W)
            du .= f.(u) .+ p * W
        end
        t0 = 0.0
        tf = 5.0
        tspan = (t0, tf)
        p = 0.01
        n = 5
        u0 = 0.1 .+ 0.7rand(n)
        prob = RODEProblem(f!, u0, tspan, p, rand_prototype=zeros(n))
        tspanbwd = (tf, t0)
        probbwd = remake(prob, tspan = (tf, t0))
        solfwd = @test_nowarn solve(prob, RandomHeun(), dt=1/100)
        solbwd = @test_nowarn solve(probbwd, RandomHeun(), dt=1/100)

        M = function (du, u, p, t, W)
            norm(du)
        end
        augprob = @test_nowarn augmentprob(prob, M)

        augsol = @test_nowarn solve(augprob, RandomHeun(), dt=1/100, save_noise=true)
        @test length(augsol.u[begin]) == 2n + 2
        @test augsol.u[begin].fwd == augsol.u[begin].bwd == u0
        @test augsol.u[begin].lfwd == augsol.u[begin].lbwd == 0.0
        @test all(≥(0.9), augsol.u[end].fwd)
        @test all(≤(0.1), augsol.u[end].bwd)
        # each direction has its own noise: the first `n` entries forward, the last `n` backward
        Wfwd = [W[1:n] for W in augsol.W.u]
        Wbwd = [W[n+1:2n] for W in augsol.W.u]
        @test augsol.u[end].lfwd ≈ trapezoid(M.(f.(getindex.(augsol.u, :fwd)) .+ p .* Wfwd, getindex.(augsol.u, :fwd), p, augsol.t, Wfwd), augsol.t) rtol = 1e-3
        @test augsol.u[end].lbwd ≈ trapezoid(M.(f.(getindex.(augsol.u, :bwd)) .+ p .* Wbwd, getindex.(augsol.u, :bwd), p, augsol.t, Wbwd), augsol.t) rtol = 1e-3
    end

    @testset "Default noise" begin
        # Without `noise` or `rand_prototype`, the noise should have the shape it has
        # when solving the original problem, and start at `t0`.
        M = function (du, u, p, t, W)
            sum(abs2, du)
        end
        t0 = 1.0
        tf = 3.0

        # scalar state: `f` and `M` receive scalar noise, stored as a one-element vector
        f = function (u, p, t, W)
            W isa Number || error("expected scalar noise, got $(typeof(W))")
            -u + W
        end
        Mscalar = function (du, u, p, t, W)
            W isa Number || error("expected scalar noise, got $(typeof(W))")
            abs2(du)
        end
        prob = RODEProblem(f, 0.5, (t0, tf))
        for direction in (:both, :forward, :backward)
            augsol = @test_nowarn solve(augmentprob(prob, Mscalar; direction), RandomHeun(), dt = 1 / 100)
            @test first(augsol.W.t) == t0
            # a copy of the noise for each direction
            @test length(augsol.W.u[end]) == (direction == :both ? 2 : 1)
        end
        lagsol = solve(LagrangianDescriptorProblem(prob, Mscalar, range(-1, 1, length = 5)), RandomHeun(), dt = 1 / 100)
        @test all(isfinite, lagsol())

        # vector state: noise with the same shape as the state
        n = 3
        f! = function (du, u, p, t, W)
            du .= -u .+ W
        end
        prob = RODEProblem(f!, rand(n), (t0, tf))
        sol = solve(prob, RandomHeun(), dt = 1 / 100)
        for direction in (:both, :forward, :backward)
            augsol = @test_nowarn solve(augmentprob(prob, M; direction), RandomHeun(), dt = 1 / 100)
            @test first(augsol.W.t) == t0
            @test size(sol.W.u[end]) == (n,)
            # a copy of the noise for each direction
            @test size(augsol.W.u[end]) == (direction == :both ? (2n,) : (n,))
        end
    end

    @testset "Independent backward noise" begin
        # u' = W: the forward solution is the integral of the forward noise, and the backward
        # solution is minus the integral of the backward noise
        M = function (du, u, p, t, W)
            abs(du)
        end
        endpoints(prob, trajectories) = solve(
            EnsembleProblem(augmentprob(prob, M); output_func = (sol, i) -> (sol.u[end], false)),
            RandomHeun(), EnsembleSerial(); trajectories, dt = 1 / 100,
        ).u

        # default noise: independent noise for each direction, so the endpoints are uncorrelated
        prob = RODEProblem((u, p, t, W) -> W, 0.0, (0.0, 1.0))
        ends = endpoints(prob, 2000)
        @test abs(cor([u.fwd for u in ends], [u.bwd for u in ends])) < 0.1

        # a given noise process: both directions use it, so the backward solution is minus
        # the forward one, with the same descriptor
        prob = RODEProblem((u, p, t, W) -> W, 0.0, (0.0, 1.0); noise = WienerProcess(0.0, 0.0))
        ends = endpoints(prob, 20)
        @test all(u -> u.bwd ≈ -u.fwd && u.lbwd ≈ u.lfwd, ends)
    end
end

@testset "Lag Desc RODE" begin
    @testset "linear OOP RODEProblem" begin
        f = function (u, p, t, W)
            u + p * W
        end
        p = 0.1
        t0 = 0.0
        tf = 5.0
        tspan = (t0, tf)
        u0 = 0.5
        prob = RODEProblem(f, u0, tspan, p)

        M = function (du, u, p, t, W)
            sum(abs2, du)
        end
        uu0 = range(-1.0, 1.0, length = 101)
        lagprob = @test_nowarn LagrangianDescriptorProblem(prob, M, uu0)
        lagsol = @test_nowarn solve(lagprob, RandomHeun(), dt=1/100)
        @test_nowarn solve(lagprob, RandomHeun(), dt=1/100, EnsembleSerial())

        @test lagsol() == lagsol(:total) == lagsol(:both) ≈ getindex.(lagsol.enssol.u, :lfwd) + getindex.(lagsol.enssol.u, :lbwd)
        @test lagsol(:forward) ≈ getindex.(lagsol.enssol.u, :lfwd)
        @test lagsol(:backward) ≈ getindex.(lagsol.enssol.u, :lbwd)
        @test lagsol(:difference) ≈ lagsol(:forward) - lagsol(:backward)

        @test_throws ArgumentError lagsol(:blah)

        lagprob =
            @test_nowarn LagrangianDescriptorProblem(prob, M, uu0; direction = :forward)
        lagsol = @test_nowarn solve(lagprob, RandomHeun(), dt=1/100)
        @test_nowarn solve(lagprob, RandomHeun(), dt=1/100, EnsembleSerial())

        @test lagsol() ≈ getindex.(lagsol.enssol.u, :lfwd) ≈ lagsol(:forward)
        @test lagsol(:forward) ≈ getindex.(lagsol.enssol.u, :lfwd)

        lagprob =
            @test_nowarn LagrangianDescriptorProblem(prob, M, uu0; direction = :backward)
        lagsol = @test_nowarn solve(lagprob, RandomHeun(), dt=1/100)
        @test_nowarn solve(lagprob, RandomHeun(), dt=1/100, EnsembleSerial())

        @test lagsol() ≈ getindex.(lagsol.enssol.u, :lbwd) ≈ lagsol(:backward)
        @test lagsol(:backward) ≈ getindex.(lagsol.enssol.u, :lbwd)

        @test_throws ArgumentError LagrangianDescriptorProblem(
            prob,
            M,
            uu0;
            direction = :blah,
        )
    end
end
