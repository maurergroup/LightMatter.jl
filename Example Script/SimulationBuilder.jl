using Lightmatter


dim = define_sim_dimensions(Dimension=0) #For 1D call something like define_sim_dimensions(Dimension=1,Lengths=400,spacing=2.0)
las = define_laser_system(:Gaussian,FWHM=50,Power=200,hv=1.55,Transport="Optical")
sim = define_simulation_settings(nlelecphon=true,nlelecheat=true,nlphonheat=true,elecphonint=true,electemp=true,phonontemp=true) #For a Two-Temperature Model Simulation
mp = define_material_parameters(las,sim,dim,extcof=16.3,debye=162.3,noatoms=59.0,plasma=13.2,thermalcond=0.0019,
                                dos="DOS/Au_DOS.dat",secmomspecfun=1.66e-5,τf=37.7)

initialtemps=Dict("Tel"=>300.0,"Tph"=>300.0)
tspan=(-100.0, 100.0)
sys = function_builder(sim,las,dim)
sol=run_simulation(sys,initialtemps,tspan,sim,mp,las,dim,save=2.0,tolerance=1e-5,max_step=0.1,min_step=0.01)
post_production(sol,"Test.hdf5",initialtemps,"minimum",sim,mp,las,dim)