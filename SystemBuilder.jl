
function build_system(sim,mp,laser,las,con,initaltemps=Dict("Nil"=>0.0)::Dict)
    sys = get_systems(mp,sim,laser)
    connections = generate_variableconnections(sys)
    defaults = generate_parameterconnections(sys)
    u0 = generate_initalconditions(sys,initialtemps)
    p = generate_parametervalues(sys,mp,las,cons,sim)
    return connected_sys,u0
end

function equation_builder(sim,mp,laser)
    @named Phonon_temp=t_phonon_factory(mp,sim)
    @named Electron_temp = t_electron_factory(mp,sim)
    no_part = get_thermalparticles(0.0,1e-16,mp.DOS,8.617e-5)
    connections=[Electron_temp.Tph ~ Phonon_temp.Tph,
                Electron_temp.Tel ~ Phonon_temp.Tel]
    chempot = (t>=0.0) => (update_chempot!,[Electron_temp.dTel.Tel=>:Tel],
    [Electron_temp.μ=>:μ,Electron_temp.kB=>:kB],[Electron_temp.μ],(mp.DOS,no_part))
    connected = compose(ODESystem(connections,t,name=:connected,defaults=Pair{Num,Any}[Phonon_temp.kB => Electron_temp.kB
    ,Phonon_temp.λ=>Electron_temp.λ,Phonon_temp.hbar=>Electron_temp.hbar,Phonon_temp.μ=>Electron_temp.μ],discrete_events=chempot)
    ,Electron_temp,Phonon_temp)
    connected_simp=structural_simplify(connected)
    return connected_simp,Electron_temp,Phonon_temp
end

function get_systems(mp,sim,laser)
    sys = Vector{ODESystem}(undef,0)
    if sim.Systems.ElectronTemperature == true
        @named Electron_temp = t_electron_factory(mp,sim,laser)
        push!(sys,Electron_temp)
    end
    if sim.Systems.PhononTemperature == true
        @named Phonon_temp=t_phonon_factory(mp,sim)
        push!(sys,Phonon_temp)
    end
    if sim.Systems.NonEqElectrons == true

    end
    return sys
end 

function generate_variableconnections(sys)
    connections = Vector{Equation}(undef,0)
    for i in 1:length(sys)-1
        for j in i+1:length(sys)
            for x in 1:length(sys[i].unknowns)
                for y in 1:length(sys[j].unknowns)
                    if Symbol(sys[i].unknowns[x])==Symbol(sys[j].unknowns[y])
                        symstring = String(Symbol(sys[i].unknowns[x]))
                        choppedsym = Symbol(chop(symstring,tail=3))
                        push!(connections,getproperty(sys[i],choppedsym) ~ getproperty(sys[j],choppedsym))
                    end
                end
            end
        end
    end
    return connections
end

function generate_parameterconnections(sys)
    defaults = Vector{Pair{Num,Any}}(undef,0)
    for i in 1:length(sys)-1
        for j in i+1:length(sys)
            for x in 1:length(parameters(sys[i]))
                for y in 1:length(parameters(sys[j]))
                    if Symbol(parameters(sys[i])[x])==Symbol(parameters(sys[j])[y])
                        sym = Symbol(parameters(sys[i])[x])
                        push!(defaults,getproperty(sys[i],sym) => getproperty(sys[j],sym))
                    end
                end
            end
        end
    end
    return defaults
end

function generate_initalconditions(sys,initialtemps)
    u0=Vector{Pair{Num,Any}}(undef,0)
    for i in eachindex(sys)
        if nameof(sys[i]) == :Electron_temp
            push!(u0,getproperty(sys[i],:Tel)=>initialtemps["Tel"])
        end
        if nameof(sys[i]) == :Phonon_temp
            push!(u0,getproperty(sys[i],:Tph)=>initialtemps["Tph"])
        end
    end
    return u0
end

function generate_parametervalues(sys,mp,las,cons,sim)
    p=Vector{Pair{Num,Any}}(undef,0)
    for i in eachindex(sys)
        for x in parameters(sys[i])
            if !isassigned(first.p,x)
                
    return p
end

function run_dynamics(connected_eq,Tel_eq,Tph_eq,las,mp)
    u0=[Tel_eq.Tel=>300.0,
        Tph_eq.Tph=>300.0]
    p=[Tel_eq.kB => 8.617e-5,
    Tel_eq.μ => 0.0,
    Tel_eq.λ => mp.λ,
    Tel_eq.hbar => 0.6582,
    Tel_eq.FWHM => las.FWHM,
    Tel_eq.ϵ => mp.ϵ,
    Tel_eq.ϕ => las.Power,
    Tel_eq.R => las.Reflectivity,
    Tph_eq.n => mp.n,
    Tph_eq.θ => mp.θ]

    tspan=(-250.0,250.0)
    prob=ODEProblem(connected_eq,u0,tspan,p)
    sol=solve(prob,Tsit5();abstol=1e-3,reltol=1e-3)
    return sol
end