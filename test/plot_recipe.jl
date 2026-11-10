using LagrangianDescriptors
using OrdinaryDiffEq, StochasticDiffEq
using Plots, Test
using SciMLBase: EnsembleSerial, remake
using LagrangianDescriptors: augmentprob
using LinearAlgebra: norm

@testset "Plot recipes" begin
    @testset "Duffing" begin
        f! = function (du, u, p, t)
            x, y = u
            A, ω = p
            du[1] = y
            du[2] = x - x^3 + A * cos(ω * t)
        end
        
        u0 = [0.5, 2.2]
        tspan = (0.0, 13.0)
        A = 0.3; ω = π; p = (A, ω)
        
        prob = ODEProblem(f!, u0, tspan, p)

        M = function (du, u, p, t) norm(du) end

        uu0 = [[x, y] for y in range(-1.0, 1.0, length=301), x in range(-1.8, 1.8, length=301)]

        lagprob = @test_nowarn LagrangianDescriptorProblem(prob, M, uu0)

        lagsol = @test_nowarn solve(lagprob, Tsit5())

        @test_nowarn plot(lagsol)
        @test_nowarn plot(lagsol, :both)
        @test_nowarn plot(lagsol, :forward)
        @test_nowarn plot(lagsol, :backward)
        @test_nowarn plot(lagsol, :difference)
        @test_throws ArgumentError plot(lagsol, :blah)

        # unsupported collection of initial conditions: a vector of two-element vectors
        lagsol = solve(LagrangianDescriptorProblem(prob, M, [[x, 0.0] for x in range(-1.8, 1.8, length=11)]), Tsit5())
        @test_throws ArgumentError plot(lagsol)
    end

    @testset "Cubic" begin
        f = function (u, p, t)
            u - u^3
        end
        t0 = 0.0
        tf = 5.0
        tspan = (t0, tf)
        u0 = 0.5
        prob = ODEProblem(f, u0, tspan)

        M = function (du, u, p, t)
            norm(du)
        end
        
        uu0 = range(-1.0, 1.0, length=101)

        lagprob = @test_nowarn LagrangianDescriptorProblem(prob, M, uu0)

        lagsol = @test_nowarn solve(lagprob, Tsit5())

        @test_nowarn plot(lagsol)
        @test_nowarn plot(lagsol, :both)
        @test_nowarn plot(lagsol, :forward)
        @test_nowarn plot(lagsol, :backward)
        @test_nowarn plot(lagsol, :difference)
        @test_throws ArgumentError plot(lagsol, :blah)
    end

    @testset "Realizations" begin
        # averaged solutions are plotted like any other
        f = function (u, p, t)
            [u[2], u[1] - u[1]^3]
        end
        g = function (u, p, t)
            [0.0, 0.1]
        end
        prob = SDEProblem(f, g, [0.0, 0.0], (0.0, 2.0))
        uu0 = [[x, y] for y in range(-1, 1, length = 11), x in range(-1, 1, length = 11)]
        lagsol = solve(LagrangianDescriptorProblem(prob, IncrementDescriptor(Δt = 0.05), uu0), EM(), dt = 0.05; realizations = 3)
        @test length(lagsol()) == length(uu0)
        @test_nowarn plot(lagsol)
        @test_nowarn plot(lagsol, :forward)
    end
end
