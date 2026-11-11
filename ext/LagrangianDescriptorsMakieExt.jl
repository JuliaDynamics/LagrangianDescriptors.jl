module LagrangianDescriptorsMakieExt

using LagrangianDescriptors: LagrangianDescriptorSolution, _plot_data
using Makie: Makie

Makie.used_attributes(::Type{<:Makie.Plot}, ::LagrangianDescriptorSolution, args...) = (:idxs,)

function Makie.plottype(lagsol::LagrangianDescriptorSolution, args...)
    kind = _plot_data(lagsol, lagsol.direction).kind
    return kind == :path ? Makie.Lines : kind == :heatmap ? Makie.Heatmap : Makie.Scatter
end

function Makie.convert_arguments(
        PT::Type{<:Makie.Plot}, lagsol::LagrangianDescriptorSolution,
        direction::Symbol = lagsol.direction; idxs = (1, 2),
    )
    data = _plot_data(lagsol, direction, idxs)
    if data.kind == :path
        return Makie.convert_arguments(PT, data.x, data.z)
    elseif data.kind == :heatmap
        # Makie's `heatmap(x, y, z)` has `z[i, j]` at `(x[i], y[j])`
        return Makie.convert_arguments(PT, data.x, data.y, permutedims(data.z))
    else
        return [Makie.PlotSpec(Makie.plotsym(PT), Makie.Point2f.(data.x, data.y); color = data.z)]
    end
end

end
