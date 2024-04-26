module E2TM0D
"""
Density of States of the chosen material. It takes in the DOS from a provided location and interpolates between
all points to form a function that can have any energy inserted in the form DOS(u) to generate the number of 
states at the provided energy. It assumes there are 0 states outside of the plotted DOS which is an okay 
assumption due to mathematical properties of the equation e.g. the answer will already be 0. There will be 
future development to allow for multiple DOS at different heights into the simulation.
DOS is multiplied by atomic density in nm^3 to convert from singlet atom DOS to states eV^-1 nm^-3
"""
function DensityOfStates(File)
    TotalDOS::Matrix{Float64}=readdlm(File,skipstart=4)
    DOS=Interpolations.interpolate(TotalDOS[:,1].+FE, TotalDOS[:,2]*n, SteffenMonotonicInterpolation())
    return DOS,TotalDOS[1,1]+FE,TotalDOS[end,1]+FE
end
"""
The Fermi Dirac Distribution at the current electronic temperature(Tel) and chemcial potential(μ).
It is used for updating the chemical potential throughout the trajectory.
"""
function FermiDirac(E::Float64,μ::Float64,Tel)::Float64
    return 1/(exp((E-μ)/(kB*Tel))+1)
end
"""
Calculates the number of electrons in the system at the given chemical potential(μ).
This should be equal to the number of electrons capture by the DOS e.g. for Cu it is 11.
This can then be checked against in the ChemicalPotential function to update the chemical 
potential, μ for the current electronic temperature Tel.
"""
function NumberOfElectrons(μ::Float64,Tel::Float64,DOS,lb,ub)::Float64
    int(u,p)=DOS(u)*FermiDirac(u,μ,Tel)
    prob=IntegralProblem(int,lb,ub)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10)
    return sol.u
end
"""
Calculates the chemical potential(μ). The function currently calculates the chemcial potential
at the given temperature by finding a root in the equation No.elec-currentno.elec=0 where the 
No.elec is the correct value and the current is calculated from NumberOfElectrons. 
"""
function ChemicalPotential(μ,Tel,DOS,DOSne,lb,ub)
    x0=μ
    f(x)=DOSne-NumberOfElectrons(x,Tel,DOS,lb,ub)
    prob=ZeroProblem(f,x0)
    return solve(prob,Order16();atol=1e-10,rtol=1e-10)
end
"""
Calculates the nonlinear electronic heat capacity for the system. This depends on the current electronic temperature(Tel)
like the linear case but also depends on the denisty of states and the rate of change of the Fermi-Dirac 
distribution with respect to temperature. The final units are eV K^-1 nm^-3
"""
function dFDdt(E,μ,Tel)
    A=(E-μ)/kB
    Frac=exp(A/Tel)/Tel^2
    Fermi=(exp(A/Tel)+1)^-2
    if isnan(A*Frac*Fermi)
        return 0
    else
        return A*Frac*Fermi
    end
end

function Cel_nonlinear(μ,Tel,DOS,lb,ub)
    int(u,p)=ForwardDiff.derivative(x->(exp((u-μ)/(kB*x))+1)^-1,Tel)*DOS(u)*u
    prob=IntegralProblem(int,lb,ub)
    sol=solve(prob,HCubatureJL(initdiv=3);reltol=1e-5,abstol=1e-5)
    return sol.u::Float64
end
"""
This is a crude approximation that can lead to temperatures as far as 25% below what they should be according
to the non linear approximation. Therefore it is recommended to perform the calculation with the non-linear form
as it is less parameterised and is not much slower than the linear form of the electronic heat capacity calculated
here. The units are eV K^-1 nm^-3
"""
function Cel_linear(Tel) 
    cel =γ*Tel #linear heat capacity of electorns according to Sommerfeld model
    return cel
end
"""
The heat capacity of the lattice is calculated here. This is a non-linear form with minimal approximations 
using composite Simpson's rule similar to that of Oscar's 2TM. The final units are eV K^-1 nm^-3
"""
function Cph(Tph) 
    uplim=θ./Tph #Generates an upperlimit for each z position in the Tph
    int(u,p)=(u^4*exp(u))/((exp(u)-1)^2)
    prob = IntegralProblem(int,0.0,uplim)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10) #sovles the integral within Cph for each position in Tph
    cph::Float64=(9*n*kB)*(Tph/θ)^3 *sol.u # Finds Cph for each z coordinate by broadcasting equation over Tph and integral
    return cph
end 

"""
What this equals is never really given a term but it a prefactor that links the electron-electorn lifetime to the
plasma frequency. Again it can be precomputed before the dynamics. This would then capture the height dependence
if implemented into the plasma frequency. The final unit is fs
"""
function tau_0()
    τ=128/(sqrt(3)*pi^2*ω)
    return τ
end
"""
Non-equilibrium electron-thermalised electron lifetime which shows the palsma frequency dependence through τ
It also has an energy which peaks at the Fermi level resulting in electrons at the Fermi Level having a much longer
lifetime than those around it. The final unit is still fs
"""
function tau_ee(τ,μ,Tel,E)
    τee=τ*μ^2 /((E-μ)^2 +(pi*kB*Tel)^2)
    return τee
end
"""
The non-equilibrium electron lifetime due to electorn-phonon scattering. The inverse is the rate of electron-phonon
scattering. This requires knowledge of the quasiparticle free-flight time through the system so requires DFT 
calcuations to acquire one of the parameters. This can be approximated to relative accuracy with τth which
is the lifetime of thermalised quasi-particles. The unit is fs
"""
function tau_ep()
    τep=τf*hv/(kB*θ)
    return τep
end
"""
An approximation to the electron-phonon lifetime used when the quasiparticle free-flight time isn't known. 
This is a rather parameterised approximation that has an error of around 30% which is not fantastic. However 
the alternative requires expensive DFT calculations to perform whereas the specific ehat capacity and 
electron-phonon coupling constants of materials are well known. The unit is fs.
"""
function tau_th(g,Tel,Tph)
    τth=(γ*(Tel+Tph)/(2*g))
    return τth
end
"""
The partial derivative of the Fermi function with respect to energy. This is an analytical expression so can 
be solved exactly. This is utilised within the more accurate temeprature-dependent electron-phonon coupling 
constant. The isnan solution is a problem as at low temps the electron-phonon coupling ≠ 0.
"""
function dFDdE(μ,E,Tel)
    v=1/(kB*Tel)
    u=-exp((E-μ)/(kB*Tel))
    w=(exp.((E-μ)/(kB*Tel))+1)^2
    f=v*u/w
    if isnan(f)
        f=0.0
    end
    return f
end
"""
A temperature-dependent electron-phonon coupling constant. The value is much different to that of the 
temp-independent version that can be utilised. This does require knowledge of some additional parameters
that can be found in https://0-journals-aps-org.pugwash.lib.warwick.ac.uk/prb/pdf/10.1103/PhysRevB.77.075133.
The final units are eV fs^-1 K^-1 nm^-3
"""
function gep_nonlinear(Tel,μ,DOS)
    int(u,p)=DOS(u)^2*-dFDdE(μ,u,Tel)
    prob=IntegralProblem(int,-20,10)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-5,abstol=1e-5)
    λ=29e-6
    gep::Float64=sol.u*(pi*kB*λ)/(DOS(μ)*hbar)
    return gep
end
"""
Transfer of heat from the electronic bath to the phonon bath. g can be found either from a constant or the non-linear
form calculated in gep_nonlinear. This results in an equation with the units eV fs^-1 nm^-3
"""
function GTelTph(g,Tel,Tph,μ,Input,DOS)
    if Input["Gep"]=="L"
            G=g*(Tel-Tph)
    elseif Input["Gep"]=="NL"
        g=gep_nonlinear(Tel,μ,DOS)
        G=g*(Tel-Tph)
    end
    return G
end
"""
The strength of the laser calculated from the FUll-Width at Half-Maximum (s) as well as the fluence(ϕ)
and th time (tprime). The laser pulse is offset by 100fs as seen in the (tprime-100) component of the 
exponential. The result is a pulse strenght in eV nm^-2 fs^-1.
"""
function LaserStrength(tprime)
    return α*sqrt(4*log(2)/pi)*(ϕ/FWHM)*exp(-4*log(2)*((tprime-Delay)/FWHM)^2)
end
function gaussian(x)
    return 1/(σ*sqrt(2*pi))*exp(-(x-hv)^2/(2*σ^2))
end

function probint(DOS,E_h,E_e,μ,Tel)
    if (E_h-E_e)>0
        probability = DOS(E_h)*DOS(E_e)*FermiDirac(E_h,μ,Tel)*(1-FermiDirac(E_e,μ,Tel))*gaussian(E_h-E_e)
    else
        probability = DOS(E_h)*DOS(E_e)*FermiDirac(E_h,μ,Tel)*(1-FermiDirac(E_e,μ,Tel))*gaussian(E_e-E_h)
    end
    return probability
end

function particlenumber(vec,DOS,ERange)
    DOSvec=DOS.(ERange)
    dE=ERange[2]-ERange[1]
    return sum(vec.*DOSvec.*dE)
end

function particleconstant(DOS,μ,Tel)
    ERange=range(-4+μ,4+μ,step=0.01)
    elec=probint.(Ref(DOS),ERange.-hv,ERange,μ,Tel)
    hol=probint.(Ref(DOS),ERange,ERange.+hv,μ,Tel)
    f(x)=(particlenumber(elec*x,DOS,ERange))-particlenumber(hol,DOS,ERange)
    x0=1
    prob=ZeroProblem(f,x0)
    return solve(prob,Order16();atol=1e-10,rtol=1e-10)
end

function FDChange(E,DOS,μ,Tel,particle)
    elec=probint(DOS,E-hv,E,μ,Tel)
    hol=probint(DOS,E,E+hv,μ,Tel)
    return (particle*elec).-hol
end
"""
This function determines the number of non-equilibrium electrons generated at each time step by
dividing the laser pulse by an integral that depends on the Fermi function and the density of states.
The integral should ideally be over infinite energies but due to the fact that this would require a 
very large DOS and that all the electrons and holes will be around the Fermi level, a range of ±3hv
appears to be appropriate. The units are fs^-1
"""
function stepint(μ,Tel,DOS,particle)
    int(u,p)=FDChange(u,DOS,μ,Tel,particle)*DOS(u)*u
    prob=IntegralProblem(int,μ-(3*hv),μ+(3*hv))
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-7,abstol=1e-7)
    return sol.u
end

function stepsize(μ,Tel,tprime,DOS,particle)
    return LaserStrength(tprime)/stepint(μ,Tel,DOS,particle)
end
"""
Calculates the internal energy transfer from the non-equilibrium electrons to the thermalised electrons over
an infinitesmially small time. This is done via an integral over all energy of the inverse of the
non-equilibrium electron lifetime due to electorn-electorn scattering multiplied by the DOS*energy and 
then the relaxed fermi dirac step change. This is the change calculated in ActualFDChange multiplied 
by an exponential of the difference between the relaxed time and current time over both lifetimes.
The units are eV fs^-2 nm^-3
"""
function ueeintegral!(u,μ,Tel,t,tprime,DOS,τ,particle)
    uee=(((u-μ)^2+(pi*kB*Tel)^2)/(τ*μ^2))*FDChange(u,DOS,μ,Tel,particle)*exp((((u-μ)^2+(pi*kB*Tel)^2)/μ^2)*-(t-tprime)/τ)*DOS(u)*u #1/τee*FDChange*exp(-((t-tprime)/τee)-((t-tprime)/τep))*DOS(E)*E
    return uee
end

function instantuee(τ,μ,Tel,tprime,t,τep,DOS,particle)
    step=stepsize(μ,Tel,tprime,DOS,particle)
    int(u,p)=ueeintegral!(u,μ,Tel,t,tprime,DOS,τ,particle)
    prob=IntegralProblem(int,μ-(3*hv),μ+(3*hv))
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-7,abstol=1e-7)
    return (step)*sol.u
end
"""
The total transfer of energy from the non-equilibrium electrons to the electronic bath due to electron-
electron scattering. This is an integral over all previous time of the instaneous transfer. This
generates an energy transfer with units of eV fs^-1 nm^-3
"""
function totalUee(g,Tph,τ,μ,Tel,t,DOS,Input,particle)
    if Input["Gep"]=="L"
        g=g
    elseif Input["Gep"]=="NL"
        g=gep_nonlinear(Tel,μ,DOS)
    end
    τep=tau_th(g,Tel,Tph)
    int(u,p)=instantuee(τ,μ,Tel,u,t,τep,DOS,particle)
    prob=IntegralProblem(int,0.0,t)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-3,abstol=1e-3)
    return sol.u::Float64
end

"""
Calculates the internal energy transfer from the non-equilibrium electrons to the phonons over
an infinitesmially small time. This is done via an integral over all energy of the inverse of the
non-equilibrium electron lifetime due to electorn-phonon scattering multiplied by the DOS*energy and 
then the relaxed fermi dirac step change. This is the change calculated in ActualFDChange multiplied 
by an exponential of the difference between the relaxed time and current time over both lifetimes.
The units are eV fs^-2 nm^-3
"""
function uepintegral!(u::Float64,τ,μ,Tel,step,t,tprime,τep,DOS)
    τee=tau_ee(τ,μ,u,Tel)
    FD=ActualFDChange(τ,μ,u,Base.steprange_last_empty)
    uep=(1/τep)*FD*exp(-((t-tprime)/τee)-((t-tprime)/p[7]))*DOS(u)*u #1/τep*FDChange*exp(-((t-tprime)/τee)-((t-tprime)/τep))*DOS(E)*E
    if isnan(uep)
        return 0.0 
    else
        return uep
    end 
end

function instantuep(τ,μ,Tel,tprime,t,τep,DOS)
    step=stepsize(μ,Tel,tprime,DOS,particle)
    int(u,p)=uepintegral!(u,τ,μ,Tel,step,t,tprime,τep,DOS)
    prob=IntegralProblem(int,μ-(3*hv),μ+(3*hv))
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-7,abstol=1e-7)
    return sol.u::Float64
end
"""
The total transfer of energy from the non-equilibrium electrons to the phonon bath due to electron-
phonon scattering. This is an integral over all previous time of the instaneous transfer. This
generates an energy transfer with units of eV fs^-1 nm^-3
"""
function totalUep(g,Tph,τ,μ,Tel,t,DOS,Input)
    if Input["Gep"]=="L"
        g=g
    elseif Input["Gep"]=="NL"
        g=gep_nonlinear(Tel,μ,DOS)
    end
    τep=tau_th(g,Tel,Tph)
    int(u,p)=instantuep(τ,μ,Tel,u,t,τep,DOS)
    prob=IntegralProblem(int,0.0,t)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-3,abstol=1e-3)
    return sol.u::Float64
end

"""
The change in electronic temperature over time. This is the sum of the thermal conductivity, electron-
phonon coupling and transfer from non-equilibrium electrons all divided by the heat capacity of the 
electronic bath. The result is K fs^-1
"""
function dTeldt(Tel,Tph,t,g,μ,Input,DOS,τ)
    if Input["Cel"]=="L"
        Ce=Cel_linear(Tel)
    elseif Input["Cel"]=="NL"
        Ce=Cel_nonlinear(μ,Tel,DOS,lb,ub)
    end
    particle=particleconstant(DOS,μ,Tel)
    #G=GTelTph(g,Tel,Tph,μ,Input,DOS)
    U=totalUee(g,Tph,τ,μ,Tel,t,DOS,Input,particle)
    Teldt=(U)/Ce
    return Teldt
end
"""
The change in phonon temperature over time. This is the sum of the electron-phonon coupling and
transfer from non-equilibrium electrons all divided by the heat capacity of the electronic bath.
The result is K fs^-1
"""
function dTphdt(Tel,Tph,t,g,μ,Input,DOS,τ)
    G=GTelTph(g,Tel,Tph,μ,Input,DOS)
    U=totalUep(g,Tph,τ,μ,Tel,t,DOS,Input)
    Tphdt=(G+U)/Cph(Tph)
    return Tphdt
end
"""
The trajectories for the coupled ODE's using DiffEq's parameter array structure. This updates
the chemcial potential at each timestep before then calculating the change in temeprature.
"""
function Trajectory!(du,u,p,t)
    if p[5]["EF"]=="NL"
        p[1]=ChemicalPotential(p[1],u[1],p[4],p[6],p[7],p[8])
    end
    println(t)
    du[1]=dTeldt(u[1],p[9],t,p[2],p[1],p[5],p[4],p[3])
    #du[2]=dTphdt(u[1],u[2],t,p[2],p[1],p[5],p[4],p[3])
end
"""
A group of parameters for the running of the trajectory.
"""
function RungeKutta(Tel,Tph,t,dt,p)
    k1=dTeldt(Tel,Tph,t,p[2],p[1],p[5],p[4],p[3])
    #l1=dTphdt(Tel,Tph,t,p[2],p[1],p[5],p[4],p[3])
    k2=dTeldt(Tel+(dt*k1/2),Tph,t+dt/2,p[2],p[1],p[5],p[4],p[3])
    #l2=dTphdt(Tel+(dt*k1/2),Tph+(dt*l1/2),t+dt/2,p[2],p[1],p[5],p[4],p[3])
    k3=dTeldt(Tel+(dt*k2/2),Tph,t+dt/2,p[2],p[1],p[5],p[4],p[3])
    #l3=dTphdt(Tel+(dt*k2/2),Tph+(dt*l2/2),t+dt/2,p[2],p[1],p[5],p[4],p[3])
    k4=dTeldt(Tel+(dt*k3),Tph,t+dt,p[2],p[1],p[5],p[4],p[3])
    #l4=dTphdt(Tel+(dt*k3),Tph+(dt*l3),t+dt,p[2],p[1],p[5],p[4],p[3])
    Tel=Tel+(dt/6*(k1+(2*k2)+(2*k3)+k4))
    #Tph=Tph+(dt/6*(l1+(2*l2)+(2*l3)+l4))
    return Tel,Tph
end
function Run0D(Input::Dict)
    τ=tau_0()
    T0=[Input["Tel"]]
    DOS,lb,ub=DensityOfStates(Input["DOS"])
    DOSne=NumberOfElectrons(FE,0.0,DOS,lb,ub)
    tspan=(0.0,Input["SimEnd"])
    μ=ChemicalPotential(FE,T0[1],DOS,DOSne,lb,ub)
    p=[μ,Input["g"],τ,DOS,Input,DOSne,lb,ub,Input["Tph"]]
#=     prob=ODEProblem(Trajectory!,T0,tspan,p)
    println("Begin Solving")
    sol=solve(prob,AutoVern7(Rodas5(autodiff=false));reltol=1e-8,abstol=1e-8,saveat=0.1,tstops=[200,250,300,350,400])
    Temp=sol.u
    Temp=reduce(vcat,transpose.(Temp))
    Telout=Temp[:,1]
    Tphout=fill(Input["Tph"],length(Telout))
    Time=sol.t =#
    Headers=["Time" "Tel" "Tph"]
    Telout=[]
    Tphout=[]
    Time=[]
    dt=0.1
    Tel=Input["Tel"]
    Tph=Input["Tph"]
    for t in range(tspan[1],tspan[2],step=dt)
        println(t)
        push!(Telout,Tel)
        push!(Tphout,Tph)
        Temps=RungeKutta(Tel,Tph,t,dt,p)
        push!(Time,t)
        Tel=Temps[1]
        Tph=Temps[2]
    end
    output=hcat(Time,Telout,Tphout)
    println("Finished Solving")
    writedlm(Input["Output"],Headers,' ')
    io=open(Input["Output"],"a")
    for row in eachrow(output)
        println(io,string(row))
    end
    close(io)
end

end