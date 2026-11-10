using LagrangianDescriptors
using OrdinaryDiffEq
using QuadGK, Test
using SciMLBase: EnsembleSerial, remake, ReturnCode, successful_retcode
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

    @testset "Non-autonomous" begin
        # Backward trajectories start at `u0` at time `t0` and go back to `2t0 - tf`.
        # A non-autonomous system with `t0 != 0` checks the time used in the backward vector field.
        t0 = 1.0
        tf = 3.0
        tspan = (t0, tf)
        tspanbwd = (t0, 2t0 - tf)
        M = function (du, u, p, t)
            norm(du)
        end
        kw = (abstol = 1e-10, reltol = 1e-10)

        # out of place
        f = function (u, p, t)
            cos(t) - u / 2
        end
        u0 = 0.5
        prob = ODEProblem(f, u0, tspan)
        solfwd = solve(prob, Tsit5(); kw...)
        solbwd = solve(remake(prob, tspan = tspanbwd), Tsit5(); kw...)
        for direction in (:both, :backward)
            augsol = solve(augmentprob(prob, M; direction), Tsit5(); kw...)
            @test augsol.u[end].bwd ≈ solbwd.u[end] rtol = 1e-6
            @test augsol.u[end].lbwd ≈ first(
                quadgk(t -> M(f(solbwd(t), nothing, t), nothing, nothing, nothing), 2t0 - tf, t0),
            ) rtol = 1e-6
        end
        augsol = solve(augmentprob(prob, M), Tsit5(); kw...)
        @test augsol.u[end].fwd ≈ solfwd.u[end] rtol = 1e-6
        @test augsol.u[end].lfwd ≈ first(
            quadgk(t -> M(f(solfwd(t), nothing, t), nothing, nothing, nothing), t0, tf),
        ) rtol = 1e-6

        # in place
        f! = function (du, u, p, t)
            du .= cos(t) .- u ./ 2
        end
        u0 = [0.5, -0.2]
        prob = ODEProblem(f!, u0, tspan)
        solbwd = solve(remake(prob, tspan = tspanbwd), Tsit5(); kw...)
        for direction in (:both, :backward)
            augsol = solve(augmentprob(prob, M; direction), Tsit5(); kw...)
            @test augsol.u[end].bwd ≈ solbwd.u[end] rtol = 1e-6
        end
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

        @test lagsol() == lagsol(:total) == lagsol(:both) ≈ getindex.(lagsol.enssol.u, :lfwd) + getindex.(lagsol.enssol.u, :lbwd)
        @test lagsol(:forward) ≈ getindex.(lagsol.enssol.u, :lfwd)
        @test lagsol(:backward) ≈ getindex.(lagsol.enssol.u, :lbwd)
        @test lagsol(:difference) ≈ lagsol(:forward) - lagsol(:backward)

        @test_throws ArgumentError lagsol(:blah)

        lagprob =
            @test_nowarn LagrangianDescriptorProblem(prob, M, uu0; direction = :forward)
        lagsol = @test_nowarn solve(lagprob, Tsit5())
        @test_nowarn solve(lagprob, Tsit5(), EnsembleSerial())

        @test lagsol() ≈ getindex.(lagsol.enssol.u, :lfwd) ≈ lagsol(:forward)
        @test lagsol(:forward) ≈ getindex.(lagsol.enssol.u, :lfwd)

        lagprob =
            @test_nowarn LagrangianDescriptorProblem(prob, M, uu0; direction = :backward)
        lagsol = @test_nowarn solve(lagprob, Tsit5())
        @test_nowarn solve(lagprob, Tsit5(), EnsembleSerial())

        @test lagsol() ≈ getindex.(lagsol.enssol.u, :lbwd) ≈ lagsol(:backward)
        @test lagsol(:backward) ≈ getindex.(lagsol.enssol.u, :lbwd)

        @test_throws ArgumentError LagrangianDescriptorProblem(
            prob,
            M,
            uu0;
            direction = :blah,
        )
    end

    @testset "Dissipative system with max_trajectory_value" begin
        # Damped Duffing oscillator - dissipative system where backward integration diverges
        function damped_duffing!(du, u, p, t)
            x, y = u
            γ = p  # damping coefficient
            du[1] = y
            du[2] = x - x^3 - γ * y
        end

        γ = 0.3  # damping
        tspan = (0.0, 100.0)
        prob = ODEProblem(damped_duffing!, [0.5, 0.5], tspan, γ)

        uu0 = [[x, y] for y in range(-2, 2, 5), x in range(-2, 2, 5)]
        M = (du, u, p, t) -> norm(du)
        lagprob = LagrangianDescriptorProblem(prob, M, uu0; direction=:both)

        # Solve with termination threshold
        lagsol = @test_nowarn solve(lagprob, Tsit5(); max_trajectory_value=1e3)

        # Check that terminated field exists and is a BitVector
        @test lagsol.terminated isa BitVector
        @test length(lagsol.terminated) == length(uu0)

        # Check that we can still extract LD values
        ld_total = @test_nowarn lagsol(:total)
        ld_fwd = @test_nowarn lagsol(:forward)
        ld_bwd = @test_nowarn lagsol(:backward)

        @test length(ld_total) == length(uu0)
        @test all(isfinite, ld_fwd)  # Forward should be finite (dissipative attracts)

        # Some trajectories should have terminated (backward diverges for dissipative systems)
        @test any(lagsol.terminated)

        # Test without max_trajectory_value (default Inf) - terminated should all be false.
        # Only integrate forward in time, where the dissipative system stays bounded;
        # backward trajectories would blow up and run until the solver's maxiters.
        lagprob_fwd = LagrangianDescriptorProblem(prob, M, uu0; direction=:forward)
        lagsol_no_term = @test_nowarn solve(lagprob_fwd, Tsit5())
        @test lagsol_no_term.terminated isa BitVector
        @test !any(lagsol_no_term.terminated)

        # Return codes: trajectories stopped by the threshold are `Terminated`
        @test lagsol.retcodes isa Vector{ReturnCode.T}
        @test length(lagsol.retcodes) == length(uu0)
        @test lagsol.terminated == (lagsol.retcodes .== ReturnCode.Terminated)
        @test all(==(ReturnCode.Success), lagsol_no_term.retcodes)
        @test all(s -> s.retcode isa ReturnCode.T, lagsol.enssol.u)

        # Return codes: diverging backward trajectories that hit `maxiters` are flagged,
        # even though they are not `terminated`
        lagprob_bwd = LagrangianDescriptorProblem(prob, M, uu0[1, :]; direction=:backward)
        lagsol_maxiters = solve(lagprob_bwd, Tsit5(); maxiters = 1000)
        @test any(==(ReturnCode.MaxIters), lagsol_maxiters.retcodes)
        @test !all(successful_retcode, lagsol_maxiters.retcodes)
        @test !any(lagsol_maxiters.terminated)
    end
end
