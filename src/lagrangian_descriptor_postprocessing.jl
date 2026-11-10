"""
    lagrangian_descriptor(sol::ODESolution, M)

Computes the Lagrangian descriptor of a single trajectory of an `ODEProblem`.

It takes the provided infinitesimal Lagrangian descriptor of the form `M=M(du, u, p, t)`, along with the time interval given in `sol.prob.tspan`, to compute
```math
L = \\int_{t_0}^{t_f} M(du(t), u(t), p, t) \\;\\mathrm{d}t,
```
where ``t_0`` and ``t_f`` are the initial and final times in `sol.prob.tspan`, `u(t) = sol(t)`, `p = sol.prob.p`, and `du(t)` is computed from `sol` in a way dependent on whether `sol.prob` is in place or out of place.

The time integration is done via `QuadGK.quadgk`.
"""
function lagrangian_descriptor(sol::ODESolution, M)
    t0, tf = extrema(sol.prob.tspan)
    if isinplace(sol.prob)
        u = similar(sol.prob.u0)
        du = similar(sol.prob.u0)
        integrand = function (t)
            sol(u, t) # interpolate in place
            sol.prob.f(du, u, sol.prob.p, t)
            return M(du, u, sol.prob.p, t)
        end
    else
        integrand = function (t)
            u = sol(t)
            du = sol.prob.f(u, sol.prob.p, t)
            return M(du, u, sol.prob.p, t)
        end
    end
    return first(quadgk(integrand, t0, tf))
end

"""
    lagrangian_descriptor(sol::RODESolution, M)

Computes the Lagrangian descriptor of a single trajectory of a `RODEProblem` or of an `SDEProblem` (whose solutions are also `RODESolution`s).

It takes the provided infinitesimal Lagrangian descriptor in the form `M=M(du, u, p, t, W)`, along with the time interval given in `sol.prob.tspan`, to compute
```math
L = \\int_{t_0}^{t_f} M(du(t), u(t), p, t, W(t)) \\;\\mathrm{d}t,
```
where ``t_0`` and ``t_f`` are the initial and final times in `sol.prob.tspan`, `u(t) = sol(t)`, `p = sol.prob.p`, `W(t) = first(sol.W(t))`, and `du(t)` is either `du = sol.prob.f(sol(t), sol.prob.p, t, first(sol.W(t)))` or given by `sol.prob.f(du, sol(t), sol.prob.p, t, first(sol.W(t)))`, depending on whether `sol.prob` is in place or out of place.

For an `SDEProblem`, the infinitesimal Lagrangian descriptor has the form `M=M(du, u, p, t)`, where `du` is the drift `sol.prob.f` along the solution, without the noise argument, as in [`augmentprob`](@ref) for `SDEProblem`s.

`QuadGK.quadgk` is an adaptive integration method and does not perform well on solutions of random or stochastic equations, so the integration, in these case, is performed via a simple trapezoidal rule on the mesh given by `sol.t`, which may be far com accurate depending on the mesh.
"""
function lagrangian_descriptor(sol::RODESolution, M)
    # The possible dependence of M on the noise makes things difficult for the adaptive quadgk, so just use a direct trapezoidal rule.
    t0, tf = sol.prob.tspan
    res = 0.0
    du = isinplace(sol.prob) ? similar(sol.prob.u0) : copy(sol.prob.u0)
    ti1 = sol.t[1]
    mi1 = 0.0
    for i in eachindex(sol.t)
        ti = sol.t[i]
        dt = ti - ti1
        mi = _local_descriptor(sol, M, du, i)
        res += (mi + mi1) * dt / 2
        ti1 = ti
        mi1 = mi
    end
    return (tf ≥ t0 ? 1 : -1) * res
end

# Infinitesimal Lagrangian descriptor at the `i`-th time of a RODE or SDE solution: for a
# RODE, `f` and `M` also take the noise; for an SDE, `f` is the drift and takes no noise.
function _local_descriptor(sol, M, du, i)
    prob = sol.prob
    u, t = sol.u[i], sol.t[i]
    args = prob isa SDEProblem ? () : (sol.W.W[i],)
    if isinplace(prob)
        prob.f(du, u, prob.p, t, args...)
    else
        du = prob.f(u, prob.p, t, args...)
    end
    return M(du, u, prob.p, t, args...)
end
