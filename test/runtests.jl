using OrdinaryDiffEq, LagrangianDescriptors, Test
using LinearAlgebra: norm
using QuadGK: quadgk

@testset "Linear OOP ODEProblem" begin
    f = function (u, p, t) p * u end
    p = 0.1
    t0 = 0.0
    tf = 5.0
    tspan = (t0, tf)
    u0 = 0.5
    prob = ODEProblem(f, u0, tspan, p)

    M = function (du, u, p, t) sum(abs2, du) end

    augprob = @test_nowarn augmentprob(prob, M)
    augsol = @test_nowarn solve(augprob, Tsit5())
    @test length(augsol.u[begin]) == 4
    @test augsol.u[begin].fwd == augsol.u[begin].bwd == u0
    @test augsol.u[end].fwd ≈ u0 * exp(p * tf)
    @test augsol.u[end].bwd ≈ u0 * exp(-p * tf)
    @test augsol.u[begin].lfwd ≈ 0.0
    @test augsol.u[end].lfwd ≈ p * u0^2 * (exp(2 * p * tf) - 1) / 2 atol=0.01
    @test augsol.u[begin].lbwd ≈ 0.0
    @test augsol.u[end].lbwd ≈ p * u0^2 * (1 - exp(-2 * p * tf)) / 2 atol=0.01

    augprob = @test_nowarn augmentprob(prob, M; direction=:forward)
    augsol = @test_nowarn solve(augprob, Tsit5())
    @test length(augsol.u[begin]) == 2
    @test augsol.u[begin].fwd == u0
    @test augsol.u[end].fwd ≈ u0 * exp(p * tf)
    @test augsol.u[begin].lfwd ≈ 0.0
    @test augsol.u[end].lfwd ≈ p * u0^2 * (exp(2 * p * tf) - 1) / 2 atol=0.01

    augprob = @test_nowarn augmentprob(prob, M; direction=:backward)
    augsol = @test_nowarn solve(augprob, Tsit5())
    @test length(augsol.u[begin]) == 2
    @test augsol.u[begin].bwd == u0
    @test augsol.u[end].bwd ≈ u0 * exp(-p * tf)
    @test augsol.u[begin].lbwd ≈ 0.0
    @test augsol.u[end].lbwd ≈ p * u0^2 * (1 - exp(-2 * p * tf)) / 2 atol=0.01
end

@testset "Cubic OOP ODEProblem" begin
    f = function (u, p, t) u - u^3 end
    t0 = 0.0
    tf = 5.0
    tspan = (t0, tf)
    u0 = 0.5
    prob = ODEProblem(f, u0, tspan)
    tspanbwd = (tf, t0)
    probbwd = remake(prob, tspan = (tf, t0))
    solfwd = @test_nowarn solve(prob, Tsit5())
    solbwd = @test_nowarn solve(probbwd, Tsit5())

    M = function (du, u, p, t) norm(du) end
    augprob = @test_nowarn augmentprob(prob, M)

    augsol = @test_nowarn solve(augprob, Tsit5())
    @test length(augsol.u[begin]) == 4
    @test augsol.u[begin].fwd == augsol.u[begin].bwd == u0
    @test augsol.u[begin].lfwd == augsol.u[begin].lbwd == 0.0
    @test augsol.u[end].fwd ≈ 1.0 atol = 0.01
    @test augsol.u[end].bwd ≈ 0.0 atol = 0.01
    @test augsol.u[end].lfwd ≈ first(quadgk(t -> M(f(solfwd(t), nothing, t), nothing, nothing, nothing), 0.0, 5.0)) atol = 0.01
    @test augsol.u[end].lbwd ≈ first(quadgk(t -> M(f(solbwd(t), nothing, t), nothing, nothing, nothing), 0.0, 5.0)) atol = 0.01
end
