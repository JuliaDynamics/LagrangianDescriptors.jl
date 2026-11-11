using LagrangianDescriptors
using OrdinaryDiffEq, StochasticDiffEq
using QuadGK, Test
using SciMLBase: EnsembleSerial, EnsembleProblem, remake, ReturnCode
using LagrangianDescriptors: augmentprob
using LinearAlgebra: norm
using Statistics: mean, var, cor
using DiffEqNoiseProcess: WienerProcess

# endpoints of the forward and backward paths of an ensemble of augmented solutions
function endpoints(augprob, alg, trajectories; kwargs...)
    ensprob = EnsembleProblem(augprob; output_func = (sol, i) -> (sol.u[end], false))
    enssol = solve(ensprob, alg, EnsembleSerial(); trajectories, kwargs...)
    return enssol.u
end

@testset "Augmented SDE" begin
    M = function (du, u, p, t)
        norm(du)
    end
    kw = (abstol = 1e-8, reltol = 1e-8)

    @testset "Zero noise" begin
        # With zero diffusion, the augmented SDE gives the same results as the augmented ODE,
        # here for a non-autonomous drift with `t0 != 0`. Euler-Maruyama with a small step is
        # used, since adaptive SDE methods are not meant for tight deterministic tolerances.
        t0, tf = 1.0, 3.0

        # out of place, scalar
        f = function (u, p, t)
            cos(t) - u / 2
        end
        g = function (u, p, t)
            0.0
        end
        odeprob = ODEProblem(f, 0.5, (t0, tf))
        sdeprob = SDEProblem(f, g, 0.5, (t0, tf))
        for direction in (:both, :forward, :backward)
            augode = solve(augmentprob(odeprob, M; direction), Tsit5(); kw...)
            augsde = solve(augmentprob(sdeprob, M; direction), EM(), dt = 1e-3)
            @test augsde.u[begin] == augode.u[begin]
            @test augsde.u[end] ≈ augode.u[end] rtol = 1e-2
        end

        # in place, vector, diagonal and non-diagonal noise
        f! = function (du, u, p, t)
            du .= cos(t) .- u ./ 2
        end
        g! = function (du, u, p, t)
            du .= 0.0
        end
        u0 = [0.5, -0.2]
        odeprob = ODEProblem(f!, u0, (t0, tf))
        sdeprob = SDEProblem(f!, g!, u0, (t0, tf))
        sdeprob_nd = SDEProblem(f!, g!, u0, (t0, tf), noise_rate_prototype = zeros(2, 3))
        for direction in (:both, :forward, :backward)
            augode = solve(augmentprob(odeprob, M; direction), Tsit5(); kw...)
            augsde = solve(augmentprob(sdeprob, M; direction), EM(), dt = 1e-3)
            @test augsde.u[end] ≈ augode.u[end] rtol = 1e-2
            augsde_nd = solve(augmentprob(sdeprob_nd, M; direction), EM(), dt = 1e-3)
            @test augsde_nd.u[end] ≈ augode.u[end] rtol = 1e-2
        end
    end

    @testset "Time-dependent noise" begin
        # dX = σ t dW from t0 = 1 to tf = 2: the forward path has variance σ² ∫_1^2 t² dt = 7σ²/3,
        # and the backward path, over [0, 1], has variance σ² ∫_0^1 t² dt = σ²/3. The forward and
        # backward noises are independent, so their endpoints are uncorrelated.
        σ = 0.5
        t0, tf = 1.0, 2.0
        trajectories = 4000
        f! = function (du, u, p, t)
            du .= 0.0
        end

        # diagonal noise
        g! = function (du, u, p, t)
            du .= σ * t
        end
        prob = SDEProblem(f!, g!, [0.0], (t0, tf))
        ends = endpoints(augmentprob(prob, M), EM(), trajectories; dt = 1 / 100)
        xfwd = [u.fwd[1] for u in ends]
        xbwd = [u.bwd[1] for u in ends]
        @test var(xfwd) ≈ 7σ^2 / 3 rtol = 0.15
        @test var(xbwd) ≈ σ^2 / 3 rtol = 0.15
        @test abs(cor(xfwd, xbwd)) < 0.1
        @test all(u -> u.lfwd == u.lbwd == 0.0, ends) # zero drift, so `M` vanishes

        # non-diagonal noise: one channel drives both components, so X - Y stays constant
        gnd! = function (dG, u, p, t)
            dG .= σ * t
        end
        u0 = [0.3, -0.1]
        prob = SDEProblem(f!, gnd!, u0, (t0, tf), noise_rate_prototype = zeros(2, 1))
        augprob = augmentprob(prob, M)
        @test size(augprob.noise_rate_prototype) == (6, 2)
        ends = endpoints(augprob, EM(), trajectories; dt = 1 / 100)
        @test all(u -> u.fwd[1] - u.fwd[2] ≈ u0[1] - u0[2], ends)
        @test all(u -> u.bwd[1] - u.bwd[2] ≈ u0[1] - u0[2], ends)
        xfwd = [u.fwd[1] for u in ends]
        xbwd = [u.bwd[1] for u in ends]
        @test var(xfwd) ≈ 7σ^2 / 3 rtol = 0.15
        @test var(xbwd) ≈ σ^2 / 3 rtol = 0.15
        @test abs(cor(xfwd, xbwd)) < 0.1
        for direction in (:forward, :backward)
            @test size(augmentprob(prob, M; direction).noise_rate_prototype) == (3, 1)
        end
    end

    @testset "Descriptor along the path" begin
        # The descriptors are the integrals of `M` of the drift along the stochastic paths.
        f! = function (du, u, p, t)
            du[1] = u[1]
            du[2] = -u[2]
        end
        g! = function (du, u, p, t)
            du .= 0.1
        end
        prob = SDEProblem(f!, g!, [0.5, 0.5], (0.0, 1.0))
        augsol = solve(augmentprob(prob, M), EM(), dt = 1 / 1000)
        trapezoid(v, t) = sum((v[begin+1:end] .+ v[begin:end-1]) .* diff(t)) / 2
        mfwd = [M([u.fwd[1], -u.fwd[2]], u.fwd, nothing, t) for (u, t) in zip(augsol.u, augsol.t)]
        mbwd = [M([-u.bwd[1], u.bwd[2]], u.bwd, nothing, t) for (u, t) in zip(augsol.u, augsol.t)]
        @test augsol.u[end].lfwd ≈ trapezoid(mfwd, augsol.t) rtol = 0.01
        @test augsol.u[end].lbwd ≈ trapezoid(mbwd, augsol.t) rtol = 0.01
    end

    @testset "Shapes and errors" begin
        f = function (u, p, t)
            u
        end
        g = function (u, p, t)
            0.1
        end
        prob = SDEProblem(f, g, 0.5, (0.0, 1.0))
        for (direction, n) in ((:both, 4), (:forward, 2), (:backward, 2))
            augsol = @test_nowarn solve(augmentprob(prob, M; direction), SOSRI())
            @test augsol.retcode == ReturnCode.Success
            @test length(augsol.u[begin]) == n
            @test all(isfinite, augsol.u[end])
        end
        @test_throws ArgumentError augmentprob(prob, M; direction = :blah)
        f2 = function (u, p, t)
            -u
        end
        g2 = function (u, p, t)
            fill(0.1, 2, 1)
        end
        prob2 = SDEProblem(f2, g2, [0.5, 0.2], (0.0, 1.0), noise_rate_prototype = zeros(2, 1))
        augsol = @test_nowarn solve(augmentprob(prob2, M), EM(), dt = 1 / 100)
        @test length(augsol.u[begin]) == 6
        @test_throws ArgumentError augmentprob(
            SDEProblem(f, g, 0.5, (0.0, 1.0), noise = WienerProcess(0.0, 0.0)), M,
        )
    end
end

@testset "Lag Desc SDE" begin
    f! = function (du, u, p, t)
        du[1] = u[1]
        du[2] = -u[2]
    end
    g! = function (du, u, p, t)
        du .= 0.1
    end
    M = function (du, u, p, t)
        norm(du)
    end
    prob = SDEProblem(f!, g!, [0.5, 0.5], (0.0, 1.0))
    uu0 = [[x, y] for y in range(-1, 1, length = 5), x in range(-1, 1, length = 5)]

    lagprob = @test_nowarn LagrangianDescriptorProblem(prob, M, uu0)
    lagsol = @test_nowarn solve(lagprob, SOSRI())
    @test length(lagsol()) == length(uu0)
    @test all(isfinite, lagsol())
    @test all(==(ReturnCode.Success), lagsol.retcodes)
    @test lagsol() ≈ lagsol(:forward) + lagsol(:backward)
    @test lagsol(:forward) ≈ getindex.(lagsol.enssol.u, :lfwd)
    @test lagsol(:backward) ≈ getindex.(lagsol.enssol.u, :lbwd)

    for direction in (:forward, :backward)
        lagsol = @test_nowarn solve(LagrangianDescriptorProblem(prob, M, uu0; direction), SOSRI())
        @test all(isfinite, lagsol())
        @test lagsol() == lagsol(direction)
    end

    # stopping diverging trajectories
    lagsol = solve(LagrangianDescriptorProblem(remake(prob, tspan = (0.0, 10.0)), M, uu0), SOSRI(); max_trajectory_value = 10.0)
    @test any(lagsol.terminated)
    @test lagsol.terminated == lagsol.terminated_forward .| lagsol.terminated_backward
    @test all(==(ReturnCode.Success), lagsol.retcodes)

    @test_throws ArgumentError LagrangianDescriptorProblem(prob, M, uu0; direction = :blah)
end

@testset "Post-proc SDE" begin
    M = function (du, u, p, t)
        norm(du)
    end
    # with zero noise, the SDE descriptor is the ODE descriptor
    f! = function (du, u, p, t)
        du .= u .- u .^ 3
    end
    g! = function (du, u, p, t)
        du .= 0.0
    end
    u0 = [0.1, 0.5]
    tspan = (0.0, 5.0)
    odesol = solve(ODEProblem(f!, u0, tspan), Tsit5())
    sdesol = solve(SDEProblem(f!, g!, u0, tspan), EM(), dt = 1 / 1000)
    @test lagrangian_descriptor(sdesol, M) ≈ lagrangian_descriptor(odesol, M) rtol = 0.01

    # with noise, it matches the augmented descriptor along the same kind of path
    g2! = function (du, u, p, t)
        du .= 0.05
    end
    sdesol = solve(SDEProblem(f!, g2!, u0, tspan), EM(), dt = 1 / 1000)
    @test isfinite(lagrangian_descriptor(sdesol, M))
    @test lagrangian_descriptor(sdesol, M) > 0
end

@testset "Stopping SDE solutions" begin
    M = function (du, u, p, t)
        norm(du)
    end
    # dX = -X dt + σ dW decays forward and grows backward in time
    for (g, kw) in (((u, p, t) -> 0.05 .* one.(u), (;)), ((u, p, t) -> fill(0.05, 1, 1), (; noise_rate_prototype = zeros(1, 1))))
        prob = SDEProblem((u, p, t) -> -u, g, [1.0], (0.0, 5.0); kw...)
        uu0 = [[x] for x in range(0.5, 2.0, length = 7)]
        stopped = solve(LagrangianDescriptorProblem(prob, M, uu0), EM(), dt = 0.01; max_trajectory_value = 5.0, shared_noise = true)
        free = solve(LagrangianDescriptorProblem(prob, M, uu0), EM(), dt = 0.01; shared_noise = true)
        @test all(stopped.terminated_backward) && !any(stopped.terminated_forward)
        @test all(stopped(:backward) .< free(:backward))
        @test all(isfinite, stopped())
    end
end

@testset "Stopping with other descriptors and problems" begin
    # increment descriptor
    prob = ODEProblem((u, p, t) -> u, 1.0, (0.0, 5.0))
    sol = solve(LagrangianDescriptorProblem(prob, IncrementDescriptor(Δt = 0.01, p = 1), [0.5, 1.0]), Tsit5(); max_trajectory_value = 10.0)
    @test all(sol.terminated_forward) && !any(sol.terminated_backward)
    @test all(10 .- [0.5, 1.0] .<= sol(:forward) .<= 10.5 .- [0.5, 1.0])
    # RODE
    rprob = RODEProblem((u, p, t, W) -> u + 0.01W, 1.0, (0.0, 5.0))
    sol = solve(LagrangianDescriptorProblem(rprob, (du, u, p, t, W) -> abs(du), [0.5, 1.0]), RandomHeun(), dt = 0.01; max_trajectory_value = 10.0)
    @test all(sol.terminated_forward) && !any(sol.terminated_backward)
    @test all(10 .- [0.5, 1.0] .<= sol(:forward) .<= 10.5 .- [0.5, 1.0])
end

@testset "Shared noise" begin
    # dX = dW: with p = 1, the increment descriptor is the sum of |ΔW| over the grid, which is
    # the same for all initial conditions if and only if they share the noise path
    prob = SDEProblem((u, p, t) -> 0.0, (u, p, t) -> 1.0, 0.0, (0.0, 1.0))
    d = IncrementDescriptor(Δt = 0.01, p = 1)
    uu0 = [0.0, 1.0, 2.0, 3.0]
    allsame(v) = all(≈(first(v)), v)

    # independent noise by default
    sol = solve(LagrangianDescriptorProblem(prob, d, uu0), EM(), dt = 0.01)
    @test !allsame(sol(:forward))
    @test !allsame(sol(:backward))

    # shared noise: the same path for all initial conditions, different in each realization
    sol = solve(LagrangianDescriptorProblem(prob, d, uu0), EM(), dt = 0.01; shared_noise = true, realizations = 3)
    values = reshape(getindex.(sol.enssol.u, :lfwd), length(uu0), 3)
    @test all(allsame, eachcol(values))
    @test !allsame(values[1, :])
    # the forward and backward paths are independent
    @test !(getindex.(sol.enssol.u, :lfwd) ≈ getindex.(sol.enssol.u, :lbwd))

    # with a seed in the problem, the noise is reproducible, independent or shared
    seeded = remake(prob, seed = 123)
    lagprob = LagrangianDescriptorProblem(seeded, d, uu0)
    sol1 = solve(lagprob, EM(), dt = 0.01)
    sol2 = solve(lagprob, EM(), dt = 0.01)
    @test sol1() == sol2()
    @test !allsame(sol1(:forward))
    shared1 = solve(lagprob, EM(), dt = 0.01; shared_noise = true)
    shared2 = solve(lagprob, EM(), dt = 0.01; shared_noise = true)
    @test shared1() == shared2()
    @test allsame(shared1(:forward))
    # without a seed, shared noise changes from one solve to the next
    @test solve(LagrangianDescriptorProblem(prob, d, uu0), EM(), dt = 0.01; shared_noise = true)() !=
        solve(LagrangianDescriptorProblem(prob, d, uu0), EM(), dt = 0.01; shared_noise = true)()

    # RODE problems: u' = W, with the default noise or a given noise process
    MW = function (du, u, p, t, W)
        abs(du)
    end
    for rprob in (RODEProblem((u, p, t, W) -> W, 0.0, (0.0, 1.0)),
                  RODEProblem((u, p, t, W) -> W, 0.0, (0.0, 1.0); noise = WienerProcess(0.0, 0.0)))
        lagprob = LagrangianDescriptorProblem(rprob, MW, uu0)
        @test !allsame(solve(lagprob, RandomHeun(), EnsembleSerial(), dt = 0.01)(:forward))
        sol = solve(lagprob, RandomHeun(), EnsembleSerial(), dt = 0.01; shared_noise = true)
        @test allsame(sol(:forward)) && allsame(sol(:backward))
    end

    # shared noise is not defined for ODE problems
    M = function (du, u, p, t)
        abs(du)
    end
    @test_throws ArgumentError solve(LagrangianDescriptorProblem(ODEProblem((u, p, t) -> u, 1.0, (0.0, 1.0)), M, uu0), Tsit5(); shared_noise = true)
end
