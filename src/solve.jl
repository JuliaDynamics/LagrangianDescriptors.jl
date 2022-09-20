function solve(prob::LagrangianDescriptorProblem, alg, args...; kwargs...)
    # This first solve was a hack on 1.7.2 rosetta; otherwise the subsequente solve would hang
    # But it is not needed on mac native and probably not on other systems as well
    # solve(prob.ensprob.prob, alg; kwargs...)
    solve(prob.ensprob, alg, args...; trajectories=length(prob.uu0), kwargs...)
end
