"""
    Assuming phonon heat capacity ≈ constant
"""
function phonon_internal_energy(mp,gc,Tph)
    return Tph * phonon_heat_capacity(mp,gc,Tph)
end

function phonon_heat_capacity(mp,gc,Tph)
    int(u,p)=(u^4*exp(u))/((exp(u)-1)^2)
    prob = IntegralProblem(int,0.0,mp.θ/Tph)
    sol=solve(prob,HCubatureJL(initdiv=2);reltol=1e-10,abstol=1e-10)
    return 9*mp.n*gc.kB*(Tph/mp.θ)^3 *sol.u
end