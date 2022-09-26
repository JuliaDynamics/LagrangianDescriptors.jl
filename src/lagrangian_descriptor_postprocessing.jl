"""
    lagrangian_descriptor(sol::ODESolution, M)

Computes the Lagrangian descriptor of a single trajectory of an `ODEProblem`.

It takes the provided infinitesimal Lagrangian descriptor `M=M(du, u, p, t)`, along with the time interval given in `sol.prob.tspan`, to compute
```math
L = \\int_{t_0}^{t_f} M(du(t), u(t), p, t) \\;\\mathrm{d}t,
```
where ``t_0`` and ``t_f`` are the initial and final times in `sol.prob.tspan`, `u(t) = sol(t)`, `p = sol.prob.p`, and `du(t)` is computed from `sol` in a way dependent on whether `sol.prob` is in place or out of place.
"""
function lagrangian_descriptor(sol::ODESolution, M)
    t0, tf = extrema(sol.prob.tspan)
    integrand = isinplace(sol.prob) ?
        function (t, du = similar(sol.prob.u0))
            sol.prob.f(du, sol(t), sol.prob.p, t)
            M(du, sol(t), sol.prob.p, t)
        end :
        function (t)
            du = sol.prob.f(sol(t), sol.prob.p, t)
            M(du, sol(t), sol.prob.p, t)
        end
    return first(quadgk(integrand, t0, tf))
end
