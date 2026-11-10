@recipe function f(lagsol::LagrangianDescriptorSolution, direction = lagsol.direction)
    _check_direction(direction, DESCRIPTOR_DIRECTIONS)
    if lagsol.uu0 isa AbstractVector{<:Number}
        seriestype --> :path
        return lagsol.uu0, lagsol(direction)
    elseif lagsol.uu0 isa AbstractMatrix{<:Vector{<:Number}} &&
           length(first(lagsol.uu0)) == 2
        seriestype --> :heatmap
        return getindex.(lagsol.uu0[1, :], 1),
        getindex.(lagsol.uu0[:, 1], 2),
        reshape(lagsol(direction), size(lagsol.uu0))
    end
end
