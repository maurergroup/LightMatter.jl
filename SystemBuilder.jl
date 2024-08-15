
function build_system(sim,mp,las,cons,dim,initialtemps=Dict("Nil"=>0.0)::Dict)
    laser=laser_factory(las,dim)
    sys = get_systems(mp,sim,laser)
    connections = generate_variableconnections(sys)
    default_params = generate_parameterconnections(sys)
    u0 = generate_initalconditions(sys,initialtemps)
    p = generate_parametervalues(sys,mp,las,cons)
    events = generate_callbacks(sim,sys,mp,cons)
    connected = compose(ODESystem(connections,t,name=:connected,defaults=default_params,discrete_events=events),sys)
    connected_sys = structural_simplify(connected)
    return connected_sys,sys,u0,p
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
                        push!(defaults,getproperty(sys[j],sym) => getproperty(sys[i],sym))
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

function generate_parametervalues(sys,mp,las,cons)
    params = Vector{Symbol}(undef,0)
    p=Vector{Pair{Num,Any}}(undef,0)
    for i in eachindex(sys)
        for x in parameters(sys[i])
            sym=Symbol(x)
            if sym ∉ params
                push!(params,sym)
                if sym in fieldnames(typeof(mp))
                    push!(p,getproperty(sys[i],sym)=>getproperty(mp,sym))
                elseif sym in fieldnames(typeof(cons))
                    push!(p,getproperty(sys[i],sym)=>getproperty(cons,sym))
                elseif sym in fieldnames(typeof(las))
                    push!(p,getproperty(sys[i],sym)=>getproperty(las,sym))
                end
            end
        end
    end
    return p
end

function generate_callbacks(sim,sys,mp,cons)
    events=[]
    if sim.ParameterApprox.ChemicalPotential == true
        if sim.Systems.NonEqElectrons==false
            for eq in sys
                if nameof(eq)==:Electron_temp
                    no_part=get_thermalparticles(mp.μ,1e-18,mp.DOS,cons.kB)
                    chempot = (t>=0.0) => (update_chempot!,[eq.dTel.Tel=>:Tel],
                    [eq.μ=>:μ,eq.kB=>:kB],[eq.μ],(mp.DOS,no_part))
                    push!(events,chempot)
                end
            end
        end
    end
    return events
end

function run_dynamics(connected_eq,u0,tspan,p)
    prob=ODEProblem(connected_eq,u0,tspan,p)
    sol=solve(prob,Tsit5();abstol=1e-3,reltol=1e-3)
    return sol
end