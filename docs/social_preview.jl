# Draws the social preview image of the repository, also shown at the top of the README and
# of the documentation, as `docs/src/img/LagrangianDescriptors_socialpreview.png`.
#
# Run it from the root of the repository, with the docs environment:
#
#     julia --project=docs docs/social_preview.jl
#
# An optional argument gives another path for the image.

using CairoMakie, OrdinaryDiffEq
using LinearAlgebra: norm
using LagrangianDescriptors

function f!(du, u, p, t)
    x, y = u
    A, ω = p
    du[1] = y
    du[2] = x - x^3 + A * sin(ω * t)
end

u0 = [0.5, 2.2]
tspan = (0.0, 13.0)
A = 5.0; ω = 2π; p = (A, ω)
prob = ODEProblem(f!, u0, tspan, p)

M(du, u, p, t) = norm(du)
uu0 = [[x, y] for y in range(-1.5, -0.5, length = 301), x in range(-0.4, 1.6, length = 601)]
lagprob = LagrangianDescriptorProblem(prob, M, uu0)

lagsol = solve(lagprob, Tsit5())

fig = Figure(size = (1280, 640), figure_padding = 0)
ax = Axis(fig[1, 1])
hidedecorations!(ax)
hidespines!(ax)
heatmap!(ax, lagsol, :forward; colormap = :inferno)
font = Makie.texfont()
text!(ax, 0.6, -1.1; text = "LagrangianDescriptors.jl", fontsize = 94, font, color = :white, align = (:center, :center))
text!(ax, 0.6, -1.25; text = "Painting the phase portrait of random and deterministic systems",
    fontsize = 36, font, color = :white, align = (:center, :center))

output = isempty(ARGS) ? joinpath(@__DIR__, "src", "img", "LagrangianDescriptors_socialpreview.png") : ARGS[1]
save(output, fig; px_per_unit = 1) # 1280 × 640, the size recommended by GitHub
println("Saved the social preview image to ", output)
