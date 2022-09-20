struct LagrangianDescriptorsSolution{T1, T2, T3}
    enssol::T1
    uu0::T2
    direction::T3
end

function (lagsol::LagrangianDescriptorsSolution)(direction::Symbol=:both)
    if direction == :both || direction == :total
        return getindex.(lagsol.enssol.u, :lfwd) + getindex.(lagsol.enssol.u, :lbwd)
    elseif direction == :forward
        return getindex.(lagsol.enssol.u, :lfwd)
    elseif direction == :backward
        return getindex.(lagsol.ensso.u, :lbwd)
    elseif direction == :difference
        return getindex.(lagsol.enssol.u, :lfwd) - getindex.(lagsol.enssol.u, :lbwd)
    else
        throw(
            ArgumentError(
                "The given argument should be a valid direction:\n" +
                "`:forward`, `:backward`, `:both` or `:total`, or `:difference`"
            )
        )
    end
end
