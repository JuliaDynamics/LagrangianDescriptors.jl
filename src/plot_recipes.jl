# Plots recipe for a `LagrangianDescriptorSolution`; see the "Plotting" page of the docs.
# Note: RecipesBase rewrites the `return`s of a recipe, so each branch below must end in its own
# `return`, in an `if`/`elseif`/`else` chain; an early `return` inside an `if` is not honored.
@recipe function f(lagsol::LagrangianDescriptorSolution, direction = lagsol.direction; idxs = (1, 2))
    data = _plot_data(lagsol, direction, idxs)
    ldlabel = direction in (:both, :total) ? "Lagrangian descriptor" :
        "$direction Lagrangian descriptor"
    if data.kind == :path
        seriestype --> :path
        xguide --> "u₀"
        yguide --> ldlabel
        legend --> false
        return data.x, data.z
    elseif data.kind == :heatmap
        seriestype --> :heatmap
        xguide --> "u₀[$(idxs[1])]"
        yguide --> "u₀[$(idxs[2])]"
        colorbar_title --> ldlabel
        return data.x, data.y, data.z
    elseif data.kind in (:surface, :scatter3d)
        # a scatter plot also for a surface, since GR, the default backend of Plots, does not
        # color a parametrized surface by other values
        seriestype --> :scatter3d
        xguide --> "u₀[$(idxs[1])]"
        yguide --> "u₀[$(idxs[2])]"
        zguide --> "u₀[$(idxs[3])]"
        # no colorbar title, which GR draws over the colorbar ticks in three dimensions
        marker_z --> vec(data.c)
        markersize --> 2
        markerstrokewidth --> 0
        label --> ""
        colorbar --> true
        return vec(data.x), vec(data.y), vec(data.z)
    else
        seriestype --> :scatter
        xguide --> "u₀[$(idxs[1])]"
        yguide --> "u₀[$(idxs[2])]"
        colorbar_title --> ldlabel
        marker_z --> data.z
        markerstrokewidth --> 0
        label --> "" # no legend entry, but keep the colorbar (`legend = false` hides both)
        colorbar --> true
        return data.x, data.y
    end
end
