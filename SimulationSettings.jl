abstract type simulation end

@kwdef struct Interactions <:simulation
    Nonequilibrium_electrons::Bool
    ElectronElectron::Bool
    ElectronPhonon::Bool
end

@kwdef struct Nonlinear <:simulation
    ChemicalPotential::Bool
    ElectronPhonon_Coupling::Bool
    ElectronHeatCapacity::Bool
end

@kwdef struct simulation_settings <: simulation
    nonlinear::Nonlinear
    interactions::Interactions
    Output_file_name::String
    Simulation_End_Time::Float64
end

function define_simulation_settings(;nlchempot=false,nlelecphon=false,
    nlelecheat=false,noneqelec=true,elecelecint=true,elecphonint=true,
    output="./Default_file.jld2",dim=0,simendtime=1000)
    
    if noneqelec==false
        if elecelecint==true
            throw(ErrorException("Electron-Electron Interactions should only 
            be true if there are non-equilibrium electrons to interact with. 
            Set elecelecint=false"))
        end
    end

    nl=Nonlinear(ChemicalPotential=nlchempot,ElectronPhononCoupling=nlelecphon,
    ElectronHeatCapacity=nlelecheat)

    interact=Interactions(NonequilibriumElectrons=noneqelec,
    ElectronElectron=elecelecint,ElectronPhonon=elecphonint)

    sim_settings=simulation_settings(nonlinear=nl,interactions=interact,
    Dimension=dim,Output_file_name=output,Simulation_End_Time=simendtime)
    return sim_settings
end


function define_simulation_settings(dict::Dict;nlchempot=dict["ChemPot"],nlelecphon=dict["ElecPhonCoup"],
    nlelecheat=dict["ElecHeatCapac"],noneqelec=dict["NoneqElec"],elecelecint=dict["Elec-Elec"],
    elecphonint=dict["Elec-Phon"],output=dict["Output"],dim=dict["Dimension"],simendtime=dict["SimTime"])

    if noneqelec==false
        if elecelecint==true
            throw(ErrorException("Electron-Electron Interactions should only 
            be true if there are non-equilibrium electrons to interact with. 
            Set elecelecint=false"))
        end
    end

    nl=Nonlinear(ChemicalPotential=nlchempot,ElectronPhononCoupling=nlelecphon,
    ElectronHeatCapacity=nlelecheat)

    interact=Interactions(NonequilibriumElectrons=noneqelec,
    ElectronElectron=elecelecint,ElectronPhonon=elecphonint)

    sim_settings=simulation_settings(nonlinear=nl,interactions=interact,
    Dimension=dim,Output_file_name=output,Simulation_End_Time=simendtime)
    return sim_settings
end