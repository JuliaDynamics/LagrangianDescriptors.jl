using LagrangianDescriptors
using OrdinaryDiffEq, StochasticDiffEq
using Plots, Test
using SciMLBase: EnsembleSerial, remake
using LagrangianDescriptors: augmentprob, _plot_data
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

        # a vector of two-element vectors is plotted as a scatter plot
        lagsol = solve(LagrangianDescriptorProblem(prob, M, [[x, 0.0] for x in range(-1.8, 1.8, length=11)]), Tsit5())
        @test plot(lagsol).series_list[1][:seriestype] == :scatter
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

    @testset "Layouts" begin
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

        # grids in both layouts: the value plotted at (x[l], y[k]) is the descriptor of the
        # initial condition [x[l], y[k]]
        for uu0 in ([[x, y] for y in ys, x in xs], [[x, y] for x in xs, y in ys])
            lagsol = solve(LagrangianDescriptorProblem(prob, M, uu0), Tsit5())
            data = _plot_data(lagsol, :both)
            @test data.kind == :heatmap
            @test data.x == xs && data.y == ys && size(data.z) == (length(ys), length(xs))
            values = Dict(zip(vec(uu0), lagsol()))
            @test all(data.z[k, l] == values[[xs[l], ys[k]]] for k in eachindex(ys), l in eachindex(xs))
            plt = plot(lagsol)
            @test plt.series_list[1][:seriestype] == :heatmap
            @test plot(lagsol; seriestype = :contourf).series_list[1][:seriestype] == :contour
        end

        # a curve of initial conditions: scatter plot colored by the descriptor
        lagsol = solve(LagrangianDescriptorProblem(prob, M, [[x, 0.3x] for x in xs]), Tsit5())
        data = _plot_data(lagsol, :forward)
        @test data.kind == :scatter && data.z == lagsol(:forward)
        plt = plot(lagsol, :forward)
        @test plt.series_list[1][:seriestype] == :scatter
        @test plt.series_list[1][:marker_z] == lagsol(:forward)

        # a scattered matrix of initial conditions, not forming a grid: scatter plot
        lagsol = solve(LagrangianDescriptorProblem(prob, M, [[x + 0.1y, y] for y in ys, x in xs]), Tsit5())
        @test _plot_data(lagsol, :both).kind == :scatter

        # a three-dimensional problem: `idxs` chooses the components of a two-dimensional slice
        g! = function (du, u, p, t)
            du[1] = u[2]
            du[2] = -u[1]
            du[3] = -u[3]
        end
        prob3 = ODEProblem(g!, zeros(3), (0.0, 1.0))
        lagsol = solve(LagrangianDescriptorProblem(prob3, M, [[x, 0.5, z] for z in ys, x in xs]), Tsit5())
        @test _plot_data(lagsol, :both, (1, 3)).kind == :heatmap
        @test _plot_data(lagsol, :both, (1, 2)).kind == :scatter # component 2 is constant
        @test plot(lagsol; idxs = (1, 3)).series_list[1][:seriestype] == :heatmap
        @test_throws ArgumentError plot(lagsol; idxs = (1, 4))
        @test_throws ArgumentError plot(lagsol; idxs = (1, 2, 4))
        @test_throws ArgumentError plot(lagsol; idxs = (1,))

        # three components: a surface, or a scatter plot for other collections; with Plots,
        # both are drawn as a 3D scatter plot
        data = _plot_data(lagsol, :forward, (3, 1, 2))
        @test data.kind == :surface
        @test size(data.x) == size(data.c) == (length(ys), length(xs))
        @test data.x == [z for z in ys, x in xs] && data.y == [x for z in ys, x in xs]
        @test vec(data.c) ≈ lagsol(:forward)
        series = plot(lagsol; idxs = (1, 2, 3)).series_list[1]
        @test series[:seriestype] == :scatter3d
        @test series[:z] ≈ [z for z in ys, x in xs][:]
        @test series[:marker_z] ≈ lagsol()
        lagsol3 = solve(LagrangianDescriptorProblem(prob3, M, [[x, x^2, -x] for x in xs]), Tsit5())
        @test _plot_data(lagsol3, :both, (1, 2, 3)).kind == :scatter3d
        @test plot(lagsol3; idxs = (1, 2, 3)).series_list[1][:seriestype] == :scatter3d

        # one-component vectors are plotted like scalars
        lagsol = solve(LagrangianDescriptorProblem(ODEProblem((u, p, t) -> -u, [0.0], (0.0, 1.0)), M, [[x] for x in xs]), Tsit5())
        @test _plot_data(lagsol, :both).kind == :path
        @test plot(lagsol).series_list[1][:seriestype] == :path
    end
end
