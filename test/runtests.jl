using LightMatter
using Test


# Some simple tests to test all completed functionality, results are saved to the test_solution.hdf5 file
# Added tests should name the section in the file the same as the test name commented above it.

io = h5open("TestSolutions.hdf5")
using LightMatter, PreallocationTools, BenchmarkTools, HDF5
egrid = collect(-6.2:0.02:6.2)
ftot = LightMatter.FermiDirac(300.0,0.0,egrid)
DOS = LightMatter.generate_DOS("../DOS/Au_DOS.dat", 1.0)
τee = 0.5 * (0.0+9.0)^2 ./((egrid.-0.0).^2 .+ (pi*Constants.kB*300.0)^2)
fneq = zeros(length(egrid))
n = LightMatter.get_thermalparticles(0.0, 1e-16, DOS, egrid)
tmp1 = DiffCache(similar(egrid))
tmp2 = DiffCache(similar(egrid))

@testset "LightMatter" begin
    #This section should contain tests that test full-code functionality,
    #E.g. test full simulation runs
end

@testset "AntennaReactor" begin
    
end

@testset "AthermalElectrons" begin
    file = io["AthermalElectrons"]
    egrid = collect(-6.2:0.02:6.2)
    ftot = LightMatter.FermiDirac(300.0,0.0,egrid)
    DOS = LightMatter.generate_DOS("../DOS/Au_DOS.dat", 1.0)
    tmp1 = DiffCache(similar(egrid))
    tmp2 = DiffCache(similar(egrid))

    LightMatter.athemexcitation!(tmp1, tmp2, ftot, egrid, DOS, 3.1, 1.0)
    @test get_tmp(tmp2,0.0) == read(file["athemexcitation"])

    τee = 0.5 * (0.0+9.0)^2 ./((egrid.-0.0).^2 .+ (pi*Constants.kB*300.0)^2)
    fneq = zeros(length(egrid))
    n = LightMatter.get_thermalparticles(0.0, 1e-16, DOS, egrid)
    LightMatter.athem_electronelectronscattering!(tmp1, tmp2, 300.0, 0.0, egrid, fneq, DOS, n, τee) 
    @test isapprox(get_tmp(tmp1, 0.0), read(file["athemelectronelectronscattering"]), atol=1e-8, rtol=1e-8)

end

@testset "DensityMatrix" begin
    
end

@testset "ElectronicDistribution" begin
    
end

@testset "ElectronicTemperature" begin
    
end

@testset "Lasers" begin
    
end

@testset "OutputProcessing" begin
    
end

@testset "PhononicDistribution" begin
    
end

@testset "PhononicTemperature" begin
    
end

@testset "PropertyFunctions" begin
    
end

@testset "SimulationTypes" begin
    
end

@testset "SystemConstruction" begin
    
end

@testset "UnitManagement" begin
    @test LightMatter.BaseUnits == (time = 1e15, length = 1e9, mass = 6.2415e30, electric_current = 1, temperature = 1, amount = 1, luminosity = 1)
    @test LightMatter.Constants == (ħ = 0.6582119569509067, kB = 8.617333262145179e-5, me = 5.685621837291226, c = 299.792458,
                                ϵ0 = 1.4185692541856925e-9, q=-0.0001602176634)
    @test LightMatter.convert_units(u"nm/fs", Unitful.c0) == LightMatter.Constants.c
end