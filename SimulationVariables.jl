function solve_NumberOfElectrons(Temp,DOS)
    int(u,p) = DOS(u) * FermiDirac(Temp,μ,u)
    prob=IntegralProblem(int,lb,ub)
    sol=solve(prob,HCubatureJL();reltol=1e-5,abstol=1e-5)
    return sol.u
end
