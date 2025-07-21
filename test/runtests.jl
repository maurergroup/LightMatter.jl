using LightMatter
using Test


# Some simple tests to test all completed functionality, results are saved to the test_solution.hdf5 file
# Added tests should name the section in the file the same as the test name commented above it.

io = h5open("TestSolutions.hdf5")

@testset "LightMatter" begin
    #This section should contain tests that test full-code functionality,
    #E.g. test full simulation runs
end

@testset "AntennaReactor" begin
    
end

@testset "AthermalElectrons" begin
    file = io["AthermalElectrons"]
    egrid = collect(-6.2:0.02:6.2)
    ftot = FermiDirac(300.0,0.0,egrid)
    DOS = get_DOS("../DOS/Au_DOS.dat")
    @test LightMatter.athemexcitation(ftot, egrid, DOS, 3.1, 1.0) == read(file["athemexcitation"])
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