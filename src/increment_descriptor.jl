"""
    IncrementDescriptor(; Δt, p = 0.5)

Lagrangian descriptor given by the sum of the `p`-th powers of the increments of the solution
over a time grid with step `Δt`, as defined for stochastic differential equations in
Balibrea-Iniesta, Lopesino, Wiggins and Mancho, *Lagrangian Descriptors for Stochastic
Differential Equations: A Tool for Revealing the Phase Portrait of Stochastic Dynamical
Systems* (2016).

It is passed in place of the infinitesimal Lagrangian descriptor `M` when building a
[`LagrangianDescriptorProblem`](@ref):

```julia
lagprob = LagrangianDescriptorProblem(prob, IncrementDescriptor(Δt = 0.05, p = 0.5), uu0)
lagsol = solve(lagprob, EM(), dt = 0.05)
```

With `tspan = (t0, tf)` and the grid `t_j = t0 + j Δt`, `j = 0, …, N`, with `t_N = tf`, the
forward Lagrangian descriptor is

```math
L_{\\mathrm{fwd}} = \\sum_{j=0}^{N-1} \\sum_k \\left| u^k(t_{j+1}) - u^k(t_j) \\right|^p,
```

where ``k`` runs over the components of the solution, and the backward one is the same sum
along the backward solution, over the grid `t0 - j Δt`, `j = 0, …, N`. Their sum is the
descriptor of the paper, over the time window `[2t0 - tf, tf]`.

The time span `tf - t0` must be a multiple of `Δt`. The solution is computed on the grid points
exactly (they are passed to the solver as `saveat` and `tstops`), not interpolated, since
interpolation would smooth out the noise. The paper uses the Euler–Maruyama method with step
`Δt`, but any algorithm can be used.

For a stochastic equation, the descriptor depends on `Δt` and does not converge as `Δt → 0`,
since the paths of a Wiener process have infinite length. The paper uses `0 < p ≤ 1`. For
`p = 1` and a differential equation, it converges to the length of the solution as `Δt → 0`.

It also works for `ODEProblem`s and `RODEProblem`s. The mean over many realizations of the
noise, also used in the paper, is obtained with the keyword `realizations` of
[`solve`](@ref).
"""
struct IncrementDescriptor{T,P}
    Δt::T
    p::P
    function IncrementDescriptor(; Δt, p = 0.5)
        Δt > 0 || throw(ArgumentError("The time step `Δt` must be positive; got `Δt = $Δt`."))
        p > 0 || throw(ArgumentError("The power `p` must be positive; got `p = $p`."))
        return new{typeof(Δt),typeof(p)}(Δt, p)
    end
end

# The time grid with step `Δt` over `tspan`, which must be a multiple of `Δt`
function _increment_grid(tspan, Δt)
    t0, tf = tspan
    N = round(Int, abs(tf - t0) / Δt)
    (N > 0 && isapprox(N * Δt, abs(tf - t0); rtol = sqrt(eps(float(Δt))))) || throw(
        ArgumentError(
            "The time span `$tspan` must be a positive multiple of `Δt = $Δt` " *
            "for an `IncrementDescriptor`.",
        ),
    )
    return range(t0, tf, length = N + 1)
end

# Sum of the `p`-th powers of the increments of the components of `x` along the saved path
function _increment_sum(path, p)
    s = zero(float(real(eltype(first(path)))))
    for j in firstindex(path):lastindex(path)-1
        a, b = path[j], path[j+1]
        for k in eachindex(a)
            s += abs(b[k] - a[k])^p
        end
    end
    return s
end

# The infinitesimal descriptor is not used, so the descriptor components stay zero
_no_descriptor(du, u, args...) = zero(eltype(u))

function _ensemble_setup(prob, d::IncrementDescriptor, direction::Symbol)
    grid = _increment_grid(prob.tspan, d.Δt)
    augprob = remake(augmentprob(prob, _no_descriptor; direction), saveat = grid, tstops = grid)
    output_func = function (sol, i)
        lfwd = direction == :backward ? nothing : _increment_sum([u.fwd for u in sol.u], d.p)
        lbwd = direction == :forward ? nothing : _increment_sum([u.bwd for u in sol.u], d.p)
        out = direction == :both ? (lfwd = lfwd, lbwd = lbwd, retcode = sol.retcode) :
              direction == :forward ? (lfwd = lfwd, retcode = sol.retcode) :
              (lbwd = lbwd, retcode = sol.retcode)
        return (out, false)
    end
    return augprob, output_func
end
