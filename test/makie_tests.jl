using LagrangianDescriptors
using OrdinaryDiffEq
using Makie, Test
# Makie functions are qualified, since Plots, loaded by other test files, exports the same names
using LinearAlgebra: norm

# the heatmap and scatter recipes wrap the actual plot, which carries the default colormap
atomic(plt) = only(Makie.collect_atomic_plots(plt))

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
            x, y, z = only(Makie.convert_arguments(Heatmap, lagsol, :forward)).args
            values = Dict(zip(vec(uu0), lagsol(:forward)))
            @test size(z) == (length(xs), length(ys))
            @test all(z[l, k] ≈ values[[xs[l], ys[k]]] for l in eachindex(xs), k in eachindex(ys))

            fig, ax, plt = Makie.plot(lagsol)
            @test atomic(plt) isa Heatmap
            @test atomic(plt).colormap[] == :inferno
            @test_nowarn Colorbar(fig[1, 2], plt)
            @test atomic(Makie.heatmap(lagsol, :backward).plot) isa Heatmap
            @test only(Makie.contourf(lagsol).plot.plots) isa Contourf
            # the default colormap gives way to the user's
            plt = Makie.heatmap(lagsol; colormap = :viridis, colorrange = (0, 1)).plot
            @test atomic(plt).colormap[] == :viridis
            @test atomic(plt).colorrange[] == [0, 1]
        end
    end

    @testset "Scatter" begin
        lagsol = solve(LagrangianDescriptorProblem(prob, M, [[x, 0.3x] for x in xs]), Tsit5())
        fig, ax, plt = Makie.plot(lagsol)
        @test atomic(plt) isa Scatter
        @test_nowarn Colorbar(fig[1, 2], plt)
        fig, ax, plt = Makie.scatter(lagsol, :forward)
        scatterplot = atomic(plt)
        @test scatterplot isa Scatter
        @test scatterplot.colormap[] == :inferno
        @test atomic(Makie.scatter(lagsol; colormap = :viridis).plot).colormap[] == :viridis
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
        @test atomic(Makie.heatmap(lagsol; idxs = (1, 3)).plot) isa Heatmap
        # the plot type depends on `idxs`: a heatmap of the components 1 and 3
        fig, ax, plt = Makie.plot(lagsol; idxs = (1, 3))
        @test atomic(plt) isa Heatmap
        @test ax isa Axis
        @test_throws ArgumentError Makie.heatmap(lagsol; idxs = (1, 4))
        @test_throws ArgumentError Makie.plot(lagsol; idxs = (1, 2, 4))
    end

    @testset "Three dimensions" begin
        g! = function (du, u, p, t)
            du[1] = u[2]
            du[2] = -u[1]
            du[3] = -u[3]
        end
        prob3 = ODEProblem(g!, zeros(3), (0.0, 1.0))
        # a matrix of initial conditions: a surface, colored by the descriptor
        uu0 = [[x, 0.5x, z] for z in ys, x in xs]
        lagsol = solve(LagrangianDescriptorProblem(prob3, M, uu0), Tsit5())
        fig, ax, plt = Makie.plot(lagsol; idxs = (1, 2, 3))
        surf = atomic(plt)
        @test surf isa Makie.Surface
        @test ax isa Axis3
        @test surf[1][] ≈ [u0[1] for u0 in uu0] && surf[3][] ≈ [u0[3] for u0 in uu0]
        @test vec(surf.color[]) ≈ lagsol()
        @test surf.colormap[] == :inferno
        @test surf.shading[] == false
        @test_nowarn Colorbar(fig[1, 2], plt)
        @test atomic(Makie.surface(lagsol; idxs = (1, 2, 3), colormap = :viridis).plot).colormap[] == :viridis
        # other plot functions plot the points
        scatterplot = atomic(Makie.scatter(lagsol, :forward; idxs = (3, 1, 2)).plot)
        @test scatterplot isa Scatter
        @test first.(scatterplot[1][]) ≈ [u0[3] for u0 in uu0][:]
        @test scatterplot.color[] ≈ lagsol(:forward)
        # any other collection of initial conditions: a scatter plot
        lagsol = solve(LagrangianDescriptorProblem(prob3, M, [[x, x^2, -x] for x in xs]), Tsit5())
        fig, ax, plt = Makie.plot(lagsol; idxs = (1, 2, 3))
        @test atomic(plt) isa Scatter
        @test ax isa Axis3
        @test last.(atomic(plt)[1][]) ≈ -xs
        @test atomic(Makie.lines(lagsol; idxs = (1, 2, 3)).plot) isa Lines
        # several solutions on the same axis
        fig = Figure()
        ax = Axis3(fig[1, 1])
        @test_nowarn Makie.plot!(ax, lagsol; idxs = (1, 2, 3))
        @test_nowarn Makie.plot!(ax, lagsol; idxs = (2, 1, 3))
    end

    @testset "Scalar" begin
        lagsol = solve(LagrangianDescriptorProblem(ODEProblem((u, p, t) -> -u, 0.0, (0.0, 1.0)), M, collect(xs)), Tsit5())
        fig, ax, plt = Makie.plot(lagsol)
        @test atomic(plt) isa Lines
        @test first.(atomic(plt)[1][]) ≈ xs
        @test last.(atomic(plt)[1][]) ≈ lagsol()
        @test atomic(Makie.scatter(lagsol, :forward).plot) isa Scatter
    end
end
