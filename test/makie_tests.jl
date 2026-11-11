using LagrangianDescriptors
using OrdinaryDiffEq
using Makie, Test
# Makie functions are qualified, since Plots, loaded by other test files, exports the same names
using LinearAlgebra: norm

@testset "Makie" begin
    @test Base.get_extension(LagrangianDescriptors, :LagrangianDescriptorsMakieExt) !== nothing

    M = function (du, u, p, t)
        norm(du)
    end
    f! = function (du, u, p, t)
        du[1] = u[2]
        du[2] = u[1] - u[1]^3 + 0.3u[2] # not symmetric, to check the orientation
    end
    prob = ODEProblem(f!, [0.0, 0.0], (0.0, 2.0))
    xs = range(-1.5, 1.5, length = 4)
    ys = range(-1.0, 1.0, length = 3)

    @testset "Grid" begin
        # in both layouts, Makie's heatmap has z[l, k] at (x[l], y[k]), which must be the
        # descriptor of the initial condition [x[l], y[k]]
        for uu0 in ([[x, y] for y in ys, x in xs], [[x, y] for x in xs, y in ys])
            lagsol = solve(LagrangianDescriptorProblem(prob, M, uu0), Tsit5())
            @test Makie.plottype(lagsol) == Heatmap
            x, y, z = Makie.convert_arguments(Heatmap, lagsol, :forward)
            values = Dict(zip(vec(uu0), lagsol(:forward)))
            @test size(z) == (length(xs), length(ys))
            @test all(z[l, k] ≈ values[[xs[l], ys[k]]] for l in eachindex(xs), k in eachindex(ys))

            fig, ax, plt = Makie.plot(lagsol)
            @test plt isa Heatmap
            @test_nowarn Colorbar(fig[1, 2], plt)
            @test Makie.heatmap(lagsol, :backward).plot isa Heatmap
            @test Makie.contourf(lagsol).plot isa Contourf
        end
    end

    @testset "Scatter" begin
        lagsol = solve(LagrangianDescriptorProblem(prob, M, [[x, 0.3x] for x in xs]), Tsit5())
        @test Makie.plottype(lagsol) == Scatter
        fig, ax, plt = Makie.plot(lagsol)
        @test_nowarn Colorbar(fig[1, 2], plt)
        fig, ax, plt = Makie.scatter(lagsol, :forward)
        scatterplot = only(Makie.collect_atomic_plots(plt))
        @test scatterplot isa Scatter
        @test scatterplot.color[] ≈ lagsol(:forward) # Makie stores colors as Float32
        @test first.(scatterplot[1][]) ≈ xs
    end

    @testset "Higher dimension and errors" begin
        g! = function (du, u, p, t)
            du[1] = u[2]
            du[2] = -u[1]
            du[3] = -u[3]
        end
        prob3 = ODEProblem(g!, zeros(3), (0.0, 1.0))
        lagsol = solve(LagrangianDescriptorProblem(prob3, M, [[x, 0.5, z] for z in ys, x in xs]), Tsit5())
        @test Makie.heatmap(lagsol; idxs = (1, 3)).plot isa Heatmap
        @test_throws ArgumentError Makie.heatmap(lagsol; idxs = (1, 4))
    end

    @testset "Scalar" begin
        lagsol = solve(LagrangianDescriptorProblem(ODEProblem((u, p, t) -> -u, 0.0, (0.0, 1.0)), M, collect(xs)), Tsit5())
        @test Makie.plottype(lagsol) == Lines
        fig, ax, plt = Makie.plot(lagsol)
        @test plt isa Lines
        @test first.(plt[1][]) ≈ xs
        @test last.(plt[1][]) ≈ lagsol()
        @test Makie.scatter(lagsol, :forward).plot isa Scatter
    end
end
