using LagrangianDescriptors
using OrdinaryDiffEq, StochasticDiffEq
using QuadGK, Test
using LagrangianDescriptors: augmentprob
using LinearAlgebra: norm

@testset "Post-process" begin
    f = function (u, p = nothing, t = nothing)
        u - u^3
    end
    f! = function (du, u, p, t)
        du .= f.(u)
    end

    t0 = 0.0
    tf = 5.0
    tspanfwd = (t0, tf)
    tspanbwd = (tf, t0)

    M = function (du, u, p, t) 
        norm(du)
    end

    @testset "out of place" begin
        u0 = 0.1 + 0.7rand()    
        probfwd = ODEProblem(f, u0, tspanfwd) 
        probbwd = remake(probfwd, tspan = tspanbwd)
        solfwd = solve(probfwd, Tsit5())
        solbwd = solve(probbwd, Tsit5())
    
        augprob = augmentprob(probfwd, M)
        augsol = solve(augprob, Tsit5())

        @test augsol.u[end].lfwd ≈ lagrangian_descriptor(solfwd, M) rtol = 0.01
        @test augsol.u[end].lbwd ≈ lagrangian_descriptor(solbwd, M) atol = 0.01
    end

    @testset "in place" begin
        n = 5
        u0 = 0.1 .+ 0.7rand(n)    
        probfwd = ODEProblem(f!, u0, tspanfwd) 
        probbwd = remake(probfwd, tspan = tspanbwd)
        solfwd = solve(probfwd, Tsit5())
        solbwd = solve(probbwd, Tsit5())
    
        augprob = augmentprob(probfwd, M)
        augsol = solve(augprob, Tsit5())

        @test augsol.u[end].lfwd ≈ lagrangian_descriptor(solfwd, M) rtol = 0.01
        @test augsol.u[end].lbwd ≈ lagrangian_descriptor(solbwd, M) atol = 0.01
    end
end

