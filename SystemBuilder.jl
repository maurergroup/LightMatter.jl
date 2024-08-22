
function build_system(sim,mp,las,cons,dim,initialtemps=Dict("Nil"=>0.0)::Dict)
    laser=laser_factory(las,dim)
    sys = generate_systems(mp,sim,laser)
    connections = generate_variableconnections(sys)
    default_params = generate_parameterconnections(sys)
    events = generate_callbacks(sim,sys,mp,cons)
    if length(sys) > 1
        connected = compose(ODESystem(connections,t,name=:connected,defaults=default_params,discrete_events=events),sys)
        connected_sys = structural_simplify(connected)
    else
        connected_sys = structural_simplify(sys[1])
    end
    #= u0 = generate_initalconditions(connected_sys,mp,cons,initialtemps)
    p = generate_parametervalues(connected_sys,mp,las,cons,initialtemps)
    precalculators!(connected_sys,p,sim,mp) =#
    return connected_sys#,u0,p
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
    for i in 1:length(sys)-1
        for j in i+1:length(sys)
            if nameof(sys[i]) == :Electron_temp && nameof(sys[j]) == :partchange
                push!(connections,getproperty(getproperty(sys[i],:dTel),:Δn) ~ D(getproperty(sys[j],:n)))
            end
            for x in 1:length(sys[i].unknowns)
                for y in 1:length(sys[j].unknowns)
                    if Symbol(sys[i].unknowns[x])==Symbol(sys[j].unknowns[y])
                        symstring = String(Symbol(sys[i].unknowns[x]))
                        if symstring[end] == ']'
                            choppedsym = Symbol(chop(symstring,tail=11,head=1))
                        else
                            choppedsym = Symbol(chop(symstring,tail=3))
                        end
                        push!(connections,getproperty(sys[i],choppedsym) ~ getproperty(sys[j],choppedsym))
                    end
                end
            end
        end
    end
    return connections
end

function generate_parameterconnections(sys)
    defaults = Vector{Pair{Union{Num,Symbolics.Arr{Num,1}},Any}}(undef,0)
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

function generate_initalconditions(sys,mp,cons,initialtemps)
    u0=Vector{Pair{Union{Num,Symbolics.Arr{Num,1}},Any}}(undef,0)
    if length(sys.systems) == 1
        if nameof(sys.systems[1]) == :dTel
            push!(u0,getproperty(sys.systems[i],:Tel)=>initialtemps["Tel"])
        elseif nameof(sys.systems[1]) == :dTph
            push!(u0,getproperty(sys.systems[i],:Tph)=>initialtemps["Tph"])
        elseif nameof(sys.systems[1]) == :dfneq
            push!(u0,getproperty(sys.systems[1],:fneq)=>@SVector zeros(length(mp.egrid)))
        end
    else
        for i in eachindex(sys.systems)
            if nameof(sys.systems[i]) == :Electron_temp
                push!(u0,getproperty(sys.systems[i],:Tel)=>initialtemps["Tel"])
            elseif nameof(sys[i]) == :Phonon_temp
                push!(u0,getproperty(sys.systems[i],:Tph)=>initialtemps["Tph"])
            elseif nameof(sys[i]) == :neq
                push!(u0,getproperty(getproperty(sys.systems[i],:dfneq),:fneq)=>@SVector zeros(length(mp.egrid)))
            elseif nameof(sys[i]) == :partchange
                push!(u0,getproperty(sys.systems[i],:n)) =>get_thermalparticles(mp.μ,1,mp.DOS,cons.kB)
            end
        end
    end
    return u0
end

function generate_parametervalues(sys,mp,las,cons,initialtemps)
    params = Vector{Symbol}(undef,0)
    p=Vector{Pair{Union{Num,Symbolics.Arr{Num,1}},Any}}(undef,0)
    if nameof(sys.systems[1]) == :dfneq 
        push!(p,getproperty(sys,:Tel)=>initialtemps["Tel"])
    end
    for x in parameters(sys)
        sym=Symbol(x)
        if sym ∉ params
            push!(params,sym)
            if sym in fieldnames(typeof(mp))
                push!(p,getproperty(sys,sym)=>getproperty(mp,sym))
            elseif sym in fieldnames(typeof(cons))
                push!(p,getproperty(sys,sym)=>getproperty(cons,sym))
            elseif sym in fieldnames(typeof(las))
                push!(p,getproperty(sys,sym)=>getproperty(las,sym))
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
            Electron_temp = sys[1]
            partchange = sys[2]
            chempot = (t>=0.0) => (update_chempotAthEM!,[Electron_temp.dTel.Tel=>:Tel,partchange.n => :n],
            [Electron_temp.μ=>:μ,Electron_temp.kB=>:kB],[Electron_temp.μ],(mp.DOS))
            push!(events,chempot)
        end
    end
    return events
end

#= function precalculators!(connected_sys,p,sim,mp)
    if 
    end 
    return p
end =#
function pre_calculations(sim,mp)
    pre_calc = Dict{String,Float64}()
    if sim.Systems.NonEqElectrons == true
        if sim.Interactions.ElectronElectron == true
            temp_dis = get_interpolate(mp.egrid,ones(length(mp.egrid)))
            n0 = get_n0(mp.DOS,mp.μ)
            u0 = get_u0(mp.DOS,mp.μ)
            pre_calc["n0"]=n0
            pre_calc["u0"]=u0
            return pre_calc
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