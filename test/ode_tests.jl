using LagrangianDescriptors
using OrdinaryDiffEq
using QuadGK, Test
using LagrangianDescriptors: augmentprob
using LinearAlgebra: norm

@testset "Augmented ODE" begin
    @testset "Linear OOP" begin
        f = function (u, p, t)
            p * u
        end
        p = 0.1
        t0 = 0.0
        tf = 5.0
        tspan = (t0, tf)
        u0 = 0.5
        prob = ODEProblem(f, u0, tspan, p)

        M = function (du, u, p, t)
            sum(abs2, du)
        end

        augprob = @test_nowarn augmentprob(prob, M)
        augsol = @test_nowarn solve(augprob, Tsit5())
        @test length(augsol.u[begin]) == 4
        @test augsol.u[begin].fwd == augsol.u[begin].bwd == u0
        @test augsol.u[end].fwd ≈ u0 * exp(p * tf)
        @test augsol.u[end].bwd ≈ u0 * exp(-p * tf)
        @test augsol.u[begin].lfwd ≈ 0.0
        @test augsol.u[end].lfwd ≈ p * u0^2 * (exp(2 * p * tf) - 1) / 2 atol = 0.01
        @test augsol.u[begin].lbwd ≈ 0.0
        @test augsol.u[end].lbwd ≈ p * u0^2 * (1 - exp(-2 * p * tf)) / 2 atol = 0.01

        augprob = @test_nowarn augmentprob(prob, M; direction = :forward)
        augsol = @test_nowarn solve(augprob, Tsit5())
        @test length(augsol.u[begin]) == 2
        @test augsol.u[begin].fwd == u0
        @test augsol.u[end].fwd ≈ u0 * exp(p * tf)
        @test augsol.u[begin].lfwd ≈ 0.0
        @test augsol.u[end].lfwd ≈ p * u0^2 * (exp(2 * p * tf) - 1) / 2 atol = 0.01

        augprob = @test_nowarn augmentprob(prob, M; direction = :backward)
        augsol = @test_nowarn solve(augprob, Tsit5())
        @test length(augsol.u[begin]) == 2
        @test augsol.u[begin].bwd == u0
        @test augsol.u[end].bwd ≈ u0 * exp(-p * tf)
        @test augsol.u[begin].lbwd ≈ 0.0
        @test augsol.u[end].lbwd ≈ p * u0^2 * (1 - exp(-2 * p * tf)) / 2 atol = 0.01

        @test_throws ArgumentError augmentprob(prob, M; direction = :blah)
    end

    @testset "Cubic OOP" begin
        f = function (u, p, t)
            u - u^3
        end
        t0 = 0.0
        tf = 5.0
        tspan = (t0, tf)
        u0 = 0.5
        prob = ODEProblem(f, u0, tspan)
        tspanbwd = (tf, t0)
        probbwd = remake(prob, tspan = (tf, t0))
        solfwd = @test_nowarn solve(prob, Tsit5())
        solbwd = @test_nowarn solve(probbwd, Tsit5())

        M = function (du, u, p, t)
            norm(du)
        end
        augprob = @test_nowarn augmentprob(prob, M)

        augsol = @test_nowarn solve(augprob, Tsit5())
        @test length(augsol.u[begin]) == 4
        @test augsol.u[begin].fwd == augsol.u[begin].bwd == u0
        @test augsol.u[begin].lfwd == augsol.u[begin].lbwd == 0.0
        @test augsol.u[end].fwd ≈ 1.0 atol = 0.01
        @test augsol.u[end].bwd ≈ 0.0 atol = 0.01
        @test augsol.u[end].lfwd ≈ first(
            quadgk(t -> M(f(solfwd(t), nothing, t), nothing, nothing, nothing), 0.0, 5.0),
        ) atol = 0.01
        @test augsol.u[end].lbwd ≈ first(
            quadgk(t -> M(f(solbwd(t), nothing, t), nothing, nothing, nothing), 0.0, 5.0),
        ) atol = 0.01
    end

    @testset "Linear IIP" begin
        f! = function (du, u, p, t)
            du .= p .* u
        end
        n = 10
        p = rand(n)
        t0 = 0.0
        tf = 5.0
        tspan = (t0, tf)
        u0 = rand(n)
        prob = ODEProblem(f!, u0, tspan, p)

        M = function (du, u, p, t)
            sum(abs2, du)
        end

        augprob = @test_nowarn augmentprob(prob, M)
        augsol = @test_nowarn solve(augprob, Tsit5())
        @test length(augsol.u[begin]) == 2n + 2
        @test augsol.u[begin].fwd == augsol.u[begin].bwd == u0
        @test augsol.u[end].fwd ≈ u0 .* exp.(p * tf) atol = 0.01
        @test augsol.u[end].bwd ≈ u0 .* exp.(-p * tf) atol = 0.01
        @test augsol.u[begin].lfwd ≈ 0.0
        @test augsol.u[end].lfwd ≈ sum(p .* u0 .^ 2 .* (exp.(2 * p * tf) .- 1) / 2) rtol =
            0.01
        @test augsol.u[begin].lbwd ≈ 0.0
        @test augsol.u[end].lbwd ≈ sum(p .* u0 .^ 2 .* (1 .- exp.(-2 * p * tf)) / 2) rtol =
            0.01

        augprob = @test_nowarn augmentprob(prob, M; direction = :forward)
        augsol = @test_nowarn solve(augprob, Tsit5())
        @test length(augsol.u[begin]) == n + 1
        @test augsol.u[begin].fwd == u0
        @test augsol.u[end].fwd ≈ u0 .* exp.(p * tf) atol = 0.01
        @test augsol.u[begin].lfwd == 0.0
        @test augsol.u[end].lfwd ≈ sum(p .* u0 .^ 2 .* (exp.(2 * p * tf) .- 1) / 2) rtol =
            0.01

        augprob = @test_nowarn augmentprob(prob, M; direction = :backward)
        augsol = @test_nowarn solve(augprob, Tsit5())
        @test length(augsol.u[begin]) == n + 1
        @test augsol.u[begin].bwd == u0
        @test augsol.u[end].bwd ≈ u0 .* exp.(-p * tf) atol = 0.01
        @test augsol.u[begin].lbwd ≈ 0.0
        @test augsol.u[end].lbwd ≈ sum(p .* u0 .^ 2 .* (1 .- exp.(-2 * p * tf)) / 2) rtol =
            0.01
    end

    @testset "Cubic IIP" begin
        f = function (u)
            u - u^3
        end
        f! = function (du, u, p, t)
            du .= f.(u)
        end
        t0 = 0.0
        tf = 5.0
        tspan = (t0, tf)
        n = 5
        u0 = 0.1 .+ 0.7rand(n)
        prob = ODEProblem(f!, u0, tspan)
        tspanbwd = (tf, t0)
        probbwd = remake(prob, tspan = (tf, t0))
        solfwd = @test_nowarn solve(prob, Tsit5())
        solbwd = @test_nowarn solve(probbwd, Tsit5())

        M = function (du, u, p, t)
            norm(du)
        end
        augprob = @test_nowarn augmentprob(prob, M)

        augsol = @test_nowarn solve(augprob, Tsit5())
        @test length(augsol.u[begin]) == 2n + 2
        @test augsol.u[begin].fwd == augsol.u[begin].bwd == u0
        @test augsol.u[begin].lfwd == augsol.u[begin].lbwd == 0.0
        @test all(≥(0.99), augsol.u[end].fwd)
        @test all(≤(0.01), augsol.u[end].bwd)
        @test augsol.u[end].lfwd ≈
              first(quadgk(t -> M(f.(solfwd(t)), nothing, nothing, nothing), 0.0, 5.0)) rtol =
            0.01
        @test augsol.u[end].lbwd ≈
              first(quadgk(t -> M(f.(solbwd(t)), nothing, nothing, nothing), 0.0, 5.0)) atol =
            0.01
    end
end

@testset "Lag Desc ODE" begin
    @testset "linear OOP ODEProblem" begin
        f = function (u, p, t)
            u
        end
        p = 0.1
        t0 = 0.0
        tf = 5.0
        tspan = (t0, tf)
        u0 = 0.5
        prob = ODEProblem(f, u0, tspan)

        M = function (du, u, p, t)
            sum(abs2, du)
        end
        uu0 = range(-1.0, 1.0, length = 101)
        lagprob = @test_nowarn LagrangianDescriptorProblem(prob, M, uu0)
        lagsol = @test_nowarn solve(lagprob, Tsit5())
        @test_nowarn solve(lagprob, Tsit5(), EnsembleSerial())

        @test lagsol() == lagsol(:total) == lagsol(:both) ≈ sum.(lagsol.enssol.u)
        @test lagsol(:forward) ≈ getindex.(lagsol.enssol.u, :lfwd)
        @test lagsol(:backward) ≈ getindex.(lagsol.enssol.u, :lbwd)
        @test lagsol(:difference) ≈ lagsol(:forward) - lagsol(:backward)

        @test_throws ArgumentError lagsol(:blah)

        lagprob =
            @test_nowarn LagrangianDescriptorProblem(prob, M, uu0; direction = :forward)
        lagsol = @test_nowarn solve(lagprob, Tsit5())
        @test_nowarn solve(lagprob, Tsit5(), EnsembleSerial())

        @test lagsol() ≈ sum.(lagsol.enssol.u) ≈ lagsol(:forward)
        @test lagsol(:forward) ≈ getindex.(lagsol.enssol.u, :lfwd)

        lagprob =
            @test_nowarn LagrangianDescriptorProblem(prob, M, uu0; direction = :backward)
        lagsol = @test_nowarn solve(lagprob, Tsit5())
        @test_nowarn solve(lagprob, Tsit5(), EnsembleSerial())

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
