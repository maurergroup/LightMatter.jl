function Parameter_factory(sim::SimulationSettings)
    @parameters μ
    μ => define_ChemicalPotential(sim.nonlinear.ChemicalPotential,mp)
end

function define_ChemicalPotential(nl::Bool,mp::MaterialParameters)
    if nl == false
        return mp.FE
    elseif nl == true
        return ChemicalPotential()
    end
end

function ChemicalPotential()

end
