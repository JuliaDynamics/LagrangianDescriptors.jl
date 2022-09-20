@recipe plot(lagsol::LagrangianDescriptorsSolution)
    lagsol.uu0, sum.(lagsol.ensprob)
end