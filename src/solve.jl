function solve(prob::LagrangianDescriptorProblem, alg, args...; kwargs...)
    solve(prob.ensprob, alg, args..., trajectories=length(prob.uu0); kwargs...)
end
