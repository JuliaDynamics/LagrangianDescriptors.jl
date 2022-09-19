function solve(prob::LagrangianDescriptorProblem, alg; kwargs...)
    # This first solve is a hack; without it, it hangs
    solve(prob.ensprob.prob, alg; kwargs...)
    solve(prob.ensprob, alg; trajectories=length(prob.uu0), kwargs...)
end
