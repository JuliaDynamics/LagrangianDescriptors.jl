module LagrangianDescriptorsMakieExt

using LagrangianDescriptors: LagrangianDescriptorSolution, _plot_data
using Makie: Makie

# The default colormap, the same as in Plots, so that the figures look the same in both;
# `colormap = ...` overrides it
const DEFAULT_COLORMAP = :inferno

Makie.used_attributes(::Type{<:Makie.Plot}, ::LagrangianDescriptorSolution, args...) = (:idxs,)

# The plot type of `plot(lagsol)` depends on `idxs`, which `Makie.plottype` does not see, so
# it is chosen in `convert_arguments`, which `plot` calls with the generic `Plot{plot}`
_default_plotsym(kind) = kind == :path ? :Lines : kind == :heatmap ? :Heatmap :
    kind == :surface ? :Surface : :Scatter

# A 3D axis for three components `idxs`; `Makie.Plot{Makie.plot}` alone does not say which
function Makie.preferred_axis_type(plot::Makie.Plot, ::LagrangianDescriptorSolution, args...)
    idxs = haskey(plot, :idxs) ? plot.idxs[] : nothing # `nothing` if not given
    return !isnothing(idxs) && length(idxs) == 3 ? Makie.Axis3 : Makie.Axis
end

function Makie.convert_arguments(
        PT::Type{<:Makie.Plot}, lagsol::LagrangianDescriptorSolution,
        direction::Symbol = lagsol.direction; idxs = (1, 2),
    )
    data = _plot_data(lagsol, direction, idxs)
    sym = PT <: Makie.Plot{Makie.plot} ? _default_plotsym(data.kind) : Makie.plotsym(PT)
    if data.kind == :path
        return [Makie.PlotSpec(sym, data.x, data.z)]
    elseif data.kind == :heatmap
        # Makie's `heatmap(x, y, z)` has `z[i, j]` at `(x[i], y[j])`
        return [Makie.PlotSpec(sym, data.x, data.y, permutedims(data.z); colormap = DEFAULT_COLORMAP)]
    elseif data.kind == :scatter
        return [Makie.PlotSpec(sym, Makie.Point2f.(data.x, data.y); color = data.z, colormap = DEFAULT_COLORMAP)]
    elseif data.kind == :surface && sym == :Surface
        # without lighting, so that the colors are those of the colorbar
        return [Makie.PlotSpec(sym, data.x, data.y, data.z; color = data.c, colormap = DEFAULT_COLORMAP, shading = false)]
    else # a surface drawn with another plot type, e.g. `scatter`, or a scatter in 3D
        return [Makie.PlotSpec(sym, Makie.Point3f.(vec(data.x), vec(data.y), vec(data.z)); color = vec(data.c), colormap = DEFAULT_COLORMAP)]
    end
end

end
