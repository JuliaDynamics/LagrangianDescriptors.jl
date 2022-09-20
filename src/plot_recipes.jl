@recipe function f(lagsol::LagrangianDescriptorSolution; direction=:both)
    if lagsol.uu0 isa AbstractVector{<:Number}
        return lagsol.uu0, lagsol(direction)
    end
end