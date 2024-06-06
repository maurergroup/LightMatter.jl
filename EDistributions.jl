using Symbolics, SymbolicNumericIntegration,ModelingToolkit,Dierckx,DelimitedFiles
using .SymbolicsInterpolation

function get_interpolate(xvals,yvals)
    return Spline1D(xvals,yvals,bc="nearest")
end
"""
    Generates an interpolation object that represents the DOS and shifted to the
    new Fermi energy
"""
function generate_DOS(File::String,FE)
    TotalDOS::Matrix{Float64}=readdlm(File,skipstart=3)
    return get_interpolate(TotalDOS[:,1].+FE,TotalDOS[:,2])
end

function FermiDirac()
    @parameters kB Tel μ
    @variables E
    return 1/(exp((E-μ)/kB*Tel)+1)
end

function dFDdT(kB,Tel,μ,E)
    e = exp((E-μ)/kB*Tel)
    return (E-μ)*e/(kB*Tel^2*(e+1)^2)
end

function integrationtest(DOS)
    @variables E
    GlobalScope(E)
    return SymbolicNumericIntegration.integrate(dFDdT()*DOS(E)*E,(E,0,20),symbolic=false)
end

DOS=generate_DOS("DOS/Au_DOS.dat",9.9)
