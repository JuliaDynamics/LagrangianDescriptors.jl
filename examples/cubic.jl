using OrdinaryDiffEq, Plots
# using LagrangianDescriptors

# Simple cubic "reaction-diffusion" equation u' = u - u^3

f(u, p, t) = u - u^3

u0 = 0.5
tspan = (0.0, 5.0)
prob = ODEProblem(f, u0, tspan)

sol = solve(prob, Tsit5())

plot(sol)

# augmented

M(du, u, p, t) = sum(abs2, du) # local descriptor

augprob = augmentprob(prob, M; direction=:backward)

augsol = solve(augprob, Tsit5())

plot(augsol)

plot(augsol, idxs=(0,1))
plot!(augsol, idxs=(0,2))
plot(augsol, idxs=(0,3,4))