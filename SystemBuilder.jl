
function build_system(sim,mp,las,cons,dim,initialtemps=Dict("Nil"=>0.0)::Dict)
    laser=laser_factory(las,dim)
    sys = generate_systems(mp,sim,laser)
    connections = generate_variableconnections(sys)
    default_params = generate_parameterconnections(sys)
    events = generate_callbacks(sim,sys,mp,cons)
    connected = compose(ODESystem(connections,t,name=:connected,defaults=default_params,discrete_events=events),sys);
    connected_sys = structural_simplify(connected)
    u0 = generate_initalconditions(connected_sys,mp,cons,initialtemps)
    p = generate_parametervalues(connected_sys,mp,las,cons,initialtemps)
    precalculators!(connected_sys,p,sim,mp)

    return connected_sys,u0,p
end

function generate_systems(mp,sim,laser)
    sys = Vector{ODESystem}(undef,0)
    if sim.Systems.ElectronTemperature == true
        @named Electron_temp = t_electron_factory(mp,sim,laser)
        push!(sys,Electron_temp)
        if sim.Systems.NonEqElectrons == true
            @named partchange = particle_change(mp.DOS,length(mp.egrid))
            push!(sys,partchange)
        end
    end
    if sim.Systems.PhononTemperature == true
        @named Phonon_temp=t_phonon_factory(mp,sim)
        push!(sys,Phonon_temp)
    end
    if sim.Systems.NonEqElectrons == true
        egl = length(mp.egrid)
        @named neq = athem_factory(sim,mp.DOS,laser,egl)
        push!(sys,neq)
    end
    return sys    
end 

function generate_variableconnections(sys)
    connections = Vector{Equation}(undef,0)
    connected = Vector{Symbol}(undef,0)
    for i in 1:length(sys)-1
        for j in i+1:length(sys)
            if nameof(sys[i]) == :Electron_temp && nameof(sys[j]) == :neq
                push!(connections,getproperty(sys[i],:relax_dis) ~ getproperty(getproperty(sys[j],:dfneq),:neqelel))
                push!(connections,getproperty(getproperty(sys[i],:dTel),:Tel) ~ getproperty(sys[j],:Tel))
            end
            if nameof(sys[i]) == :Electron_temp && nameof(sys[j]) == :partchange
                push!(connections,getproperty(getproperty(sys[i],:dTel),:Δn) ~ D(getproperty(sys[j],:n)))
            end
            for x in 1:length(sys[i].unknowns)
                sym_x = Symbol(sys[i].unknowns[x])
                for y in 1:length(sys[j].unknowns)
                    sym_y = Symbol(sys[j].unknowns[y])
                    if sym_x==sym_y
                        symstring = String(Symbol(sym_x))
                        if symstring[end] == ']'
                            choppedsym = Symbol(chop(symstring,tail=11,head=1))
                            if Symbol(getproperty(sys[j],choppedsym)) ∉ connected
                                push!(connections,getproperty(sys[i],choppedsym) ~ getproperty(sys[j],choppedsym))
                                push!(connected,Symbol(getproperty(sys[j],choppedsym)))
                            end
                        else
                            choppedsym = Symbol(chop(symstring,tail=3))
                            if Symbol(getproperty(sys[j],choppedsym)) ∉ connected
                                push!(connections,getproperty(sys[i],choppedsym) ~ getproperty(sys[j],choppedsym))
                                push!(connected,Symbol(getproperty(sys[j],choppedsym)))
                            end
                        end
                    end
                end
            end
        end
    end
    return connections
end

function generate_parameterconnections(sys)
    defaults = Vector{Pair{Union{Num,Symbolics.Arr{Num,1}},Any}}(undef,0)
    connected = Vector{Symbol}(undef,0)
    for i in 1:length(sys)-1
        for j in i+1:length(sys)
            if nameof(sys[i]) == :Electron_temp && nameof(sys[j]) == :neq
                push!(defaults,getproperty(getproperty(sys[i],:dTel),:kB) => getproperty(sys[j],:kB))
                push!(defaults,getproperty(getproperty(sys[i],:dTel),:μ) => getproperty(sys[j],:μ))
            end
            for x in 1:length(parameters(sys[i]))
                for y in 1:length(parameters(sys[j]))
                    if Symbol(parameters(sys[i])[x])==Symbol(parameters(sys[j])[y])
                        sym = Symbol(parameters(sys[i])[x])
                        if Symbol(getproperty(sys[i],sym)) ∉ connected
                            push!(defaults,getproperty(sys[j],sym) => getproperty(sys[i],sym))
                            push!(connected,Symbol(getproperty(sys[j],sym)))
                        end
                    end
                end
            end
        end
    end
    return defaults
end

function generate_initalconditions(sys,mp,cons,initialtemps)
    u0=Vector{Pair{Union{Num,Symbolics.Arr{Num,1}},Any}}(undef,0)
    if length(sys.systems) == 1
        if nameof(sys.systems[1]) == :dTel
            push!(u0,getproperty(sys.systems[i],:Tel)=>initialtemps["Tel"])
        elseif nameof(sys.systems[1]) == :dTph
            push!(u0,getproperty(sys.systems[i],:Tph)=>initialtemps["Tph"])
        elseif nameof(sys.systems[1]) == :neq
            push!(u0,getproperty(getproperty(sys.systems[1],:dfneq),:fneq)=>zeros(length(mp.egrid)))
        end
    else
        for i in eachindex(sys.systems)
            if nameof(sys.systems[i]) == :Electron_temp
                push!(u0,getproperty(getproperty(sys.systems[i],:dTel),:Tel)=>initialtemps["Tel"])
            elseif nameof(sys.systems[i]) == :Phonon_temp
                push!(u0,getproperty(sys.systems[i],:Tph)=>initialtemps["Tph"])
            elseif nameof(sys.systems[i]) == :neq
                push!(u0,getproperty(getproperty(sys.systems[i],:dfneq),:fneq)=>@SVector zeros(length(mp.egrid)))
            elseif nameof(sys.systems[i]) == :partchange
                push!(u0,getproperty(sys.systems[i],:n) => get_thermalparticles(mp.μ,1.0,mp.DOS,cons.kB))
            end
        end
    end
    return u0
end

function generate_parametervalues(sys,mp,las,cons,initialtemps)
    params = Vector{Symbol}(undef,0)
    p=Vector{Pair{Any,Any}}(undef,0)
    if length(sys.systems) == 1
        if nameof(sys.systems[1]) == :neq
            push!(p,getproperty(sys.systems[1],:Tel)=>initialtemps["Tel"])
        end
    end
    for x in parameters(sys)
        stringsym = String(Symbol(x))
        chop_idx = findfirst("₊",stringsym)[1]
        sym=Symbol(chop(stringsym,head=chop_idx,tail=0))
        if sym ∉ params
            push!(params,sym)
            if sym in fieldnames(typeof(mp))
                push!(p,x=>getproperty(mp,sym))
            elseif sym in fieldnames(typeof(cons))
                push!(p,x=>getproperty(cons,sym))
            elseif sym in fieldnames(typeof(las))
                push!(p,x=>getproperty(las,sym))
            end
        end
    end
    return p
end

function generate_callbacks(sim,sys,mp,cons)
    events=[]
    if sim.ParameterApprox.ChemicalPotential == true
        if sim.Systems.NonEqElectrons==false
            Electron_temp = sys[1]
            no_part=get_thermalparticles(mp.μ,1e-18,mp.DOS,cons.kB)
            chempot = (t>=0.0) => (update_chempotTTM!,[Electron_temp.dTel.Tel=>:Tel],
            [Electron_temp.μ=>:μ,Electron_temp.kB=>:kB],[Electron_temp.μ],(mp.DOS,no_part))
            push!(events,chempot)
        elseif sim.Systems.NonEqElectrons==true
            if sim.Systems.ElectronTemperature==true
                Electron_temp = sys[1]
                partchange = sys[2]
                chempot = (t>=0.0) => (update_chempotAthEM!,[Electron_temp.dTel.Tel=>:Tel,partchange.n => :n],
                [Electron_temp.μ=>:μ,Electron_temp.dTel.kB=>:kB],[Electron_temp.μ],(mp.DOS))
                push!(events,chempot)
            end
        end
    end
    return events
end

function precalculators!(sys,p,sim,mp)
    #push!(p,getproperty(getproperty(getproperty(sys,:Electron_temp),:dTel),:μ)=>mp.μ)
    if sim.Systems.NonEqElectrons == true
        if sim.Interactions.ElectronElectron == true
            n0 = get_n0(mp.DOS,mp.μ)
            u0 = get_u0(mp.DOS,mp.μ)
            push!(p,getproperty(getproperty(sys,:neq),:u0)=>u0)
            push!(p,getproperty(getproperty(sys,:partchange),:n0)=>n0)
        end
    else
        return nothing
    end
end

function run_dynamics(connected_eq,u0,tspan,p)
    prob=ODEProblem(connected_eq,u0,tspan,p)
    sol=solve(prob,Tsit5();abstol=1e-5,reltol=1e-5)
    return sol
end