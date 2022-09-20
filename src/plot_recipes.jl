@recipe function f(lagsol::LagrangianDescriptorSolution, direction=lagsol.direction)
    direction in (:forward, :backward, :both, :total, :difference) || throw(
        ArgumentError(
            "The given argument should be a valid direction:\n" +
            "`:forward`, `:backward`, `:both` or `:total`, or `:difference`"
        )
    )
    if lagsol.uu0 isa AbstractVector{<:Number}
        seriestype --> :path
        return lagsol.uu0, lagsol(direction)
    elseif lagsol.uu0 isa AbstractMatrix{<:Vector{<:Number}} && length(first(lagsol.uu0)) == 2
        lsol = fill(0.0, size(lagsol.uu0)...)
        if direction == :both || direction == :total
            for i in eachindex(lagsol.enssol.u) 
                lsol[i] = lagsol.enssol.u[i].lfwd + lagsol.enssol.u[i].lbwd
            end
        elseif direction == :forward
            for i in eachindex(lagsol.enssol.u) 
                lsol[i] = lagsol.enssol.u[i].lfwd
            end
        elseif direction == :backward
            for i in eachindex(lagsol.enssol.u) 
                lsol[i] = lagsol.enssol.u[i].lbwd
            end
        elseif direction == :difference
            for i in eachindex(lagsol.enssol.u) 
                lsol[i] = lagsol.enssol.u[i].lfwd - lagsol.enssol.u[i].lbwd
            end
        end
        seriestype --> :heatmap
        return getindex.(lagsol.uu0[1,:], 1), getindex.(lagsol.uu0[:,1], 2), lsol
    end
end