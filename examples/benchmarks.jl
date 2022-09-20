
using OrdinaryDiffEq, Test
using LinearAlgebra: norm
using QuadGK: quadgk
using BenchmarkTools: @btime

include("../src/LagrangianDescriptors.jl")
using .LagrangianDescriptors
using .LagrangianDescriptors: augmentprob

@testset "Augmented vs post-processing scalar ODE" begin
    f = function (u) u - u^3 end
    f! = function (du, u, p, t) du .= f.(u) end
    t0 = 0.0
    tf = 5.0
    tspan = (t0, tf)
    n = 5
    u0 = 0.1 .+ 0.7rand(n)
    prob = ODEProblem(f!, u0, tspan)
    @info "Create forward ODE problem:"
    @btime ODEProblem($f!, $u0, $tspan)
    tspanbwd = (tf, t0)
    probbwd = remake(prob, tspan = (tf, t0))
    @info "Remake for backward ODE problem:"
    @btime remake($prob, tspan = $((tf, t0)))
    solfwd = @test_nowarn solve(prob, Tsit5())
    solbwd = @test_nowarn solve(probbwd, Tsit5())
    @info "Solve forward ODE problem:"
    @btime solve($prob, $(Tsit5()))
    @info "Solve backward ODE problem:"
    @btime solve($probbwd, $(Tsit5()))

    M = function (du, u, p, t) norm(du) end
    augprob = @test_nowarn augmentprob(prob, M)
    @info "Create augmented ODE problem:"
    @btime augmentprob($prob, $M)
    augsol = @test_nowarn solve(augprob, Tsit5())
    @info "Solve Augmented ODE problem:"
    @btime solve($augprob, $(Tsit5()))

    postproc = function (sol, f, M, tspan) first(quadgk(t -> M(f.(sol(t)), nothing, nothing, nothing), first(tspan), last(tspan))) end
    @test augsol.u[end].lfwd ≈ postproc(solfwd, f, M, tspan) rtol = 0.01
    @info "Postprocessing for forward Lagrangian descriptor:"
    @btime $postproc($solfwd, $f, $M, $tspan)
    @test augsol.u[end].lbwd ≈ postproc(solbwd, f, M, tspan) atol = 0.01
    @info "Postprocessing for backward Lagrangian descriptor:"
    @btime $postproc($solbwd, $f, $M, $tspan)
end
