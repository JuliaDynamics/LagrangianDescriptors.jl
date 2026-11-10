using LagrangianDescriptors
using StochasticDiffEq
using QuadGK, Test
using SciMLBase: EnsembleSerial, remake
using LagrangianDescriptors: augmentprob
using LinearAlgebra: norm

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
        @test augsol.u[end].lfwd ≈ sum(M.(f.(getindex.(augsol.u, :fwd), p, augsol.t, first.(augsol.W.u)), getindex.(augsol.u, :fwd), p, augsol.t, first.(augsol.W.u))[begin + 1:end] .* (augsol.t[begin + 1:end] - augsol.t[begin:end-1])) rtol = 0.01
        @test augsol.u[end].lbwd ≈ sum(M.(f.(getindex.(augsol.u, :bwd), p, augsol.t, first.(augsol.W.u)), getindex.(augsol.u, :bwd), p, augsol.t, first.(augsol.W.u))[begin + 1:end] .* (augsol.t[begin + 1:end] - augsol.t[begin:end-1])) rtol = 0.01
        @test augsol.u[begin].lfwd ≈ 0.0
        @test augsol.u[begin].lbwd ≈ 0.0

        augprob = @test_nowarn augmentprob(prob, M; direction = :forward)
        augsol = @test_nowarn solve(augprob, RandomHeun(), dt=1/100, save_noise=true)
        @test length(augsol.u[begin]) == 2
        @test augsol.u[begin].fwd == u0
        @test augsol.u[begin].lfwd ≈ 0.0
        @test augsol.u[end].lfwd ≈ sum(M.(f.(getindex.(augsol.u, :fwd), p, augsol.t, first.(augsol.W.u)), getindex.(augsol.u, :fwd), p, augsol.t, first.(augsol.W.u))[begin + 1:end] .* (augsol.t[begin + 1:end] - augsol.t[begin:end-1])) rtol = 0.01

        augprob = @test_nowarn augmentprob(prob, M; direction = :backward)
        augsol = @test_nowarn solve(augprob, RandomHeun(), dt=1/100, save_noise=true)
        @test length(augsol.u[begin]) == 2
        @test augsol.u[begin].bwd == u0
        @test augsol.u[begin].lbwd ≈ 0.0
        @test augsol.u[end].lbwd ≈ sum(M.(f.(getindex.(augsol.u, :bwd), p, augsol.t, first.(augsol.W.u)), getindex.(augsol.u, :bwd), p, augsol.t, first.(augsol.W.u))[begin + 1:end] .* (augsol.t[begin + 1:end] - augsol.t[begin:end-1])) rtol = 0.01

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
        @test augsol.u[end].lfwd ≈ sum(M.(f.(getindex.(augsol.u, :fwd)) .+ p .* augsol.W.u, getindex.(augsol.u, :fwd), p, augsol.t, augsol.W.u)[begin + 1:end] .* (augsol.t[begin + 1:end] - augsol.t[begin:end-1])) rtol = 0.01
        @test augsol.u[end].lbwd ≈ sum(M.(f.(getindex.(augsol.u, :bwd)) .+ p .* augsol.W.u, getindex.(augsol.u, :bwd), p, augsol.t, augsol.W.u)[begin + 1:end] .* (augsol.t[begin + 1:end] - augsol.t[begin:end-1])) rtol = 0.01
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
            @test length(augsol.W.u[end]) == 1
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
            @test size(augsol.W.u[end]) == size(sol.W.u[end]) == (n,)
        end
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

        @test lagsol() == lagsol(:total) == lagsol(:both) ≈ sum.(lagsol.enssol.u)
        @test lagsol(:forward) ≈ getindex.(lagsol.enssol.u, :lfwd)
        @test lagsol(:backward) ≈ getindex.(lagsol.enssol.u, :lbwd)
        @test lagsol(:difference) ≈ lagsol(:forward) - lagsol(:backward)

        @test_throws ArgumentError lagsol(:blah)

        lagprob =
            @test_nowarn LagrangianDescriptorProblem(prob, M, uu0; direction = :forward)
        lagsol = @test_nowarn solve(lagprob, RandomHeun(), dt=1/100)
        @test_nowarn solve(lagprob, RandomHeun(), dt=1/100, EnsembleSerial())

        @test lagsol() ≈ sum.(lagsol.enssol.u) ≈ lagsol(:forward)
        @test lagsol(:forward) ≈ getindex.(lagsol.enssol.u, :lfwd)

        lagprob =
            @test_nowarn LagrangianDescriptorProblem(prob, M, uu0; direction = :backward)
        lagsol = @test_nowarn solve(lagprob, RandomHeun(), dt=1/100)
        @test_nowarn solve(lagprob, RandomHeun(), dt=1/100, EnsembleSerial())

        @test lagsol() ≈ sum.(lagsol.enssol.u) ≈ lagsol(:backward)
        @test lagsol(:backward) ≈ getindex.(lagsol.enssol.u, :lbwd)

        @test_throws ArgumentError LagrangianDescriptorProblem(
            prob,
            M,
            uu0;
            direction = :blah,
        )
    end
end
