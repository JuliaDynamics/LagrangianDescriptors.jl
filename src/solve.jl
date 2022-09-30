"""
    solve(prob::LagrangianDescriptorProblem, alg, args...; kwargs...)

Solve a [`LagrangianDescriptorProblem`](@ref), which amounts to solving the associated `EnsembleProblem` in `prob.ensprob` and returning a [`LagrangianDescriptorSolution`](@ref).

You should provide the necessary `args` and the desired `kwargs` for solving the associated ensemble problem for the underlying Differential Equation problem.
"""
function solve(prob::LagrangianDescriptorProblem, alg, args...; kwargs...)
    # This first solve was a hack on 1.7.2 rosetta; otherwise the subsequente solve would hang
    # But it is not needed on mac native and probably not on other systems as well
    # solve(prob.ensprob.prob, alg; kwargs...)
    ntraj = prob.method == :postprocessed ? 2length(prob.uu0) : length(prob.uu0)
    sol = solve(prob.ensprob, alg, args...; trajectories = ntraj, kwargs...)
    return LagrangianDescriptorSolution(sol, prob.uu0, prob.direction)
end
