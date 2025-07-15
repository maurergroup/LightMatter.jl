function laser(t, z, sim::Simulation)
    return temporal_laser(t, sim, Val(sim.laser.envelope)) * spatial_laser(sim, z, Val(sim.laser.Transport))
end

function temporal_laser(t,sim, type::Val{:Gaussian})
    return sqrt(4*log(2)/pi) / sim.laser.FWHM * exp(-4*log(2)*t^2/sim.laser.FWHM^2)
end

function temporal_laser(t,sim, type::Val{:HyperbolicSecant})
    return log(1+sqrt(2)) / sim.laser.FWHM * sech(2*log(1+sqrt(2))*(t/sim.laser.FWHM))^2
end

function temporal_laser(t,sim, type::Val{:Lorentzian})
    lorent = (1+(4/(1+sqrt(2)) * (t/sim.laser.FWHM)^2))^-2
    return 4 * sqrt(sqrt(2)-1) / (pi*sim.laser.FWHM) * lorent
end

function temporal_laser(t,sim, type::Val{:Rectangular})
    return -2*sim.laser.FWHM ≤ t ≤ 2*sim.laser.FWHM ? 1/(4*sim.laser.FWHM) : 0.0
end

function spatial_laser(sim::Simulation, z, type::Val{:Ballistic})
    l = sim.structure.dimension.grid[end]
    return 1 ./ (sim.laser.δb.*(1 .-exp.(-l./sim.laser.δb))) .* exp.(-z./sim.laser.δb)
end

function spatial_laser(sim::Simulation, z, type::Val{:Optical})
    l = sim.structure.dimension.grid[end]
    return 1 ./(sim.laser.ϵ .* (1 .-exp.(-l./sim.laser.ϵ))) .* exp.(-z./sim.laser.ϵ)
end

function spatial_laser(sim::Simulation, z, type::Val{:Combined})
    l = sim.structure.dimension.grid[end]
    return 1 ./ ((sim.laser.δb.+sim.laser.ϵ) * (1 .-exp.(-l./(sim.laser.δb.+sim.laser.ϵ)))) .* exp.(-z./(sim.laser.δb.+sim.laser.ϵ))
end

function laserfield(t, sim::Simulation)
    power = get_laserfield_power(sim.laser.hv)
    E_0 = laserfield_amplitude(t, power, sim, Val(sim.laser.envelope))
    E = laserfield_electric(E_0, sim.laser.hv, t)
    B = [0.0, E ./ Constants.c, 0.0]
    return (E,B)
end

function get_laserfield_power(hv::AbstractArray)
    return hv[:,2]
end

function get_laserfield_power(hv::Float64)
    return 1.0
end

function laserfield_amplitude(t, power, sim, type::Val{:Rectangular})
    return -2*sim.laser.FWHM ≤ t ≤ 2*sim.laser.FWHM ? sqrt.(2*sim.laser.ϕ.*power ./ (Constants.c*Constants.ϵ0*4*sim.laser.FWHM*sim.laser.n)) : 0.0
end

function laserfield_amplitude(t, power, sim, type::Val{:Gaussian})
    return sqrt.(2*sim.laser.ϕ*sqrt(4*log(2)).*power./(Constants.c*Constants.ϵ0*sim.laser.FWHM*sim.laser.n*sqrt(pi))).*exp(-2*log(2)*t^2/sim.laser.FWHM^2)
end

function laserfield_electric(E0, hv::AbstractArray, t)
    return [sum(E0 .* cos.(hv[:,1] * 1/(Constants.ħ*2*pi) * t)), 0.0, 0.0]
end

function laserfield_electric(E0, hv::Float64, t)
    return [sum(E0 * cos(hv * 1/(Constants.ħ*2*pi) * t)), 0.0, 0.0]
end

function photon_energytofrequency(hv)
    return hv / Constants.ħ
end

function E_magnitude(las_field, ext_field)
    sum = las_field .+ ext_field
    return sqrt((sum[1])^2 + (sum[2])^2 + (sum[3])^2)
end