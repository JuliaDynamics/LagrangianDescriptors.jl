using LagrangianDescriptors
using OrdinaryDiffEq, StochasticDiffEq
using QuadGK, Test
using SciMLBase: ReturnCode
using LinearAlgebra: norm

@testset "Increment descriptor" begin
    kw = (abstol = 1e-10, reltol = 1e-10)

    @testset "ODE, p = 1: length of the path" begin
        # For p = 1, the sum of the increments of a monotone component is its total variation.
        f = function (u, p, t)
            u
        end
        prob = ODEProblem(f, 1.0, (0.0, 1.0))
        d = IncrementDescriptor(Δt = 0.01, p = 1)
        lagsol = solve(LagrangianDescriptorProblem(prob, d, [1.0, 2.0]), Tsit5(); kw...)
        @test lagsol(:forward) ≈ [1, 2] .* (exp(1) - 1) rtol = 1e-8
        @test lagsol(:backward) ≈ [1, 2] .* (1 - exp(-1)) rtol = 1e-8
        @test lagsol() ≈ lagsol(:forward) + lagsol(:backward)
        for direction in (:forward, :backward)
            ls = solve(LagrangianDescriptorProblem(prob, d, [1.0, 2.0]; direction), Tsit5(); kw...)
            @test ls() ≈ lagsol(direction)
        end

        # in place, vector: the sum is over the components
        f! = function (du, u, p, t)
            du[1] = u[1]
            du[2] = -u[2]
        end
        prob = ODEProblem(f!, [1.0, 2.0], (0.0, 1.0))
        lagsol = solve(LagrangianDescriptorProblem(prob, d, [[1.0, 2.0]]), Tsit5(); kw...)
        @test only(lagsol(:forward)) ≈ (exp(1) - 1) + 2 * (1 - exp(-1)) rtol = 1e-8
        @test only(lagsol(:backward)) ≈ (1 - exp(-1)) + 2 * (exp(1) - 1) rtol = 1e-8
    end

    @testset "ODE, non-autonomous, t0 != 0" begin
        # u' = cos(t): forward over [1, 3] and backward over [-1, 1], so with p = 1 the sums
        # are the total variations of sin over these intervals, up to O(Δt²) at turning points.
        f = function (u, p, t)
            cos(t)
        end
        prob = ODEProblem(f, 0.0, (1.0, 3.0))
        lagsol = solve(LagrangianDescriptorProblem(prob, IncrementDescriptor(Δt = 0.01, p = 1), [0.0]), Tsit5(); kw...)
        @test only(lagsol(:forward)) ≈ first(quadgk(t -> abs(cos(t)), 1, 3)) rtol = 1e-3
        @test only(lagsol(:backward)) ≈ 2sin(1) rtol = 1e-3
    end

    @testset "ODE, p < 1" begin
        # For a smooth path, the sum is approximately Δt^(p - 1) ∫ |u'|^p dt.
        f = function (u, p, t)
            u
        end
        prob = ODEProblem(f, 1.0, (0.0, 1.0))
        Δt, p = 0.01, 0.5
        lagsol = solve(LagrangianDescriptorProblem(prob, IncrementDescriptor(; Δt, p), [1.0]), Tsit5(); kw...)
        @test only(lagsol(:forward)) ≈ Δt^(p - 1) * first(quadgk(t -> exp(t)^p, 0, 1)) rtol = 1e-4
        @test only(lagsol(:backward)) ≈ Δt^(p - 1) * first(quadgk(t -> exp(-t)^p, 0, 1)) rtol = 1e-4
    end

    @testset "SDE, Brownian motion" begin
        # For dX = σ dW, the increments are independent N(0, σ²Δt), with E|ΔX| = σ√Δt √(2/π),
        # so with p = 1 the mean of the sum over the N = T/Δt steps is N σ √Δt √(2/π), in
        # both directions.
        σ, Δt, T = 0.3, 0.02, 1.0
        N = round(Int, T / Δt)
        prob = SDEProblem((u, p, t) -> 0.0, (u, p, t) -> σ, 0.0, (0.0, T))
        lagprob = LagrangianDescriptorProblem(prob, IncrementDescriptor(Δt = Δt, p = 1), zeros(5))
        lagsol = solve(lagprob, EM(), dt = Δt; realizations = 400)
        expected = N * σ * sqrt(Δt) * sqrt(2 / π)
        @test length(lagsol()) == 5
        @test all(≈(expected; rtol = 0.05), lagsol(:forward))
        @test all(≈(expected; rtol = 0.05), lagsol(:backward))
        @test length(lagsol.retcodes) == length(lagsol.terminated) == 400 * 5
        @test all(==(ReturnCode.Success), lagsol.retcodes)
    end

    @testset "RODE" begin
        f = function (u, p, t, W)
            -u + W
        end
        prob = RODEProblem(f, 0.5, (0.0, 1.0))
        lagsol = @test_nowarn solve(LagrangianDescriptorProblem(prob, IncrementDescriptor(Δt = 0.01), range(-1, 1, length = 5)), RandomHeun(), dt = 0.01)
        @test all(isfinite, lagsol())
        @test all(>(0), lagsol())
    end

    @testset "Errors" begin
        prob = ODEProblem((u, p, t) -> u, 1.0, (0.0, 1.0))
        @test_throws ArgumentError IncrementDescriptor(Δt = 0)
        @test_throws ArgumentError IncrementDescriptor(Δt = 0.01, p = -1)
        @test_throws ArgumentError LagrangianDescriptorProblem(prob, IncrementDescriptor(Δt = 0.3), [1.0])
        @test_throws ArgumentError LagrangianDescriptorProblem(prob, IncrementDescriptor(Δt = 0.01), [1.0]; direction = :blah)
    end
end

@testset "Realizations" begin
    M = function (du, u, p, t)
        norm(du)
    end
    # deterministic problem: every realization gives the same values
    prob = ODEProblem((u, p, t) -> u - u^3, 0.5, (0.0, 2.0))
    uu0 = range(-1, 1, length = 7)
    lagprob = LagrangianDescriptorProblem(prob, M, uu0)
    lagsol1 = solve(lagprob, Tsit5())
    lagsol3 = solve(lagprob, Tsit5(); realizations = 3)
    @test length(lagsol3.enssol.u) == length(lagsol3.retcodes) == 3 * length(uu0)
    for direction in (:both, :forward, :backward, :difference)
        @test lagsol3(direction) ≈ lagsol1(direction)
    end
    @test_throws ArgumentError solve(lagprob, Tsit5(); realizations = 0)

    # stochastic problem with the integral descriptor: the mean of independent realizations
    prob = SDEProblem((u, p, t) -> -u, (u, p, t) -> 0.2, 0.5, (0.0, 1.0))
    lagprob = LagrangianDescriptorProblem(prob, M, uu0)
    lagsol = solve(lagprob, EM(), dt = 0.01; realizations = 4)
    values = reshape(getindex.(lagsol.enssol.u, :lfwd) + getindex.(lagsol.enssol.u, :lbwd), length(uu0), 4)
    @test lagsol() ≈ vec(sum(values, dims = 2)) / 4
    @test !(values[:, 1] ≈ values[:, 2]) # independent noise in each realization
end
