import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from nexusformat.nexus import *
from textwrap import wrap
from moviepy.editor import *
from moviepy.video.io.bindings import mplfig_to_npimage

mpl.rcParams['text.usetex'] = True
mpl.rcParams["mathtext.default"] = 'regular'
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
mpl.rcParams['axes.linewidth'] = 3.0

def FermiDirac(egrid,cp,T):
    k = 8.617e-5
    return 1/(np.exp((egrid-cp)/(k*T))+1)

AthEM=nxload('C:/Users/u5522838/Documents/Collaboration Docs/AthEM Nexus Files/fig6+7+8/fig6+7+8/athem_au.nxs')
BTE = nxload('C:/Users/u5522838/Documents/Collaboration Docs/AthEM Nexus Files/fig6+7+8/fig6+7+8/Boltzmann_au.nxs')

BTE_time = BTE.dynamics.time.time_discretization.nxvalue
AthEM_time = AthEM.dynamics.time.time_discretization.nxvalue

BTE_FE = BTE.materials.material_0.electrons.band_0.total.offset_energy_static.nxvalue
AthEM_FE = AthEM.materials.material_0.electrons.band_1.total.offset_energy_static.nxvalue

BTE_data=BTE.dynamics.material_0.position_0_0_0.electrons.band_0.total.distribution.distribution.nxvalue
AthEM_data=AthEM.postprocessing.total_distribution.distribution.nxvalue

AthEM_DOS = AthEM.materials.material_0.electrons.band_0.total.DOS.dos.nxvalue
BTE_DOS = BTE.materials.material_0.electrons.band_0.total.DOS.dos.nxvalue

BTE_energy = BTE.dynamics.material_0.position_0_0_0.electrons.band_0.total.distribution.energy.nxvalue+BTE_FE
AthEM_energy = AthEM.dynamics.material_0.position_0_0_0.electrons.band_0.total.distribution.energy.nxvalue+AthEM_FE


def make_frame(t):
    # sim_t = t*15*1e-15-25e-15
    # BTE_idx = min(range(len(BTE_time)), key=lambda i: abs(BTE_time[i]-sim_t))
    # AthEM_idx = min(range(len(AthEM_time)), key=lambda i: abs(AthEM_time[i]-sim_t))
    # BTE_dis=BTE_data[BTE_idx,:]
    # AthEM_dis=AthEM_data[AthEM_idx,:]
    # ax.clear()
    # ax.text(-4,1.1,str(round(sim_t/1e-15,2))+' fs',fontsize=20)
    # ax.plot(BTE_energy[200:-100],BTE_dis[200:-100],label='Boltzmann',linewidth=5,color='red')
    # ax.plot(AthEM_energy[200:-100],AthEM_dis[200:-100],label='AthEM',linewidth=5,color='blue',dashes=[6,2])
    # ax.set_xlabel(r'$\mathrm{E-E_F / eV}$',fontsize=24)
    # ax.set_ylabel('Distribution',fontsize=24)
    # ax.set_title("\n".join(wrap('Simulation with electron-electron scattering',40)),fontsize=28,fontweight='bold')
    # ax.tick_params(axis='both',which='major',labelsize=20,width=3)
    # ax.legend(prop={'size': 24})
    
    sim_t = t*10*1e-15-65e-15
    BTE_idx = min(range(len(BTE_time)), key=lambda i: abs(BTE_time[i]-sim_t))
    AthEM_idx = min(range(len(AthEM_time)), key=lambda i: abs(AthEM_time[i]-sim_t))

    BTE_cp = BTE.dynamics.material_0.position_0_0_0.electrons.band_0.total.abs_chemical_potential.abs_chemical_potential.nxvalue[BTE_idx]
    BTE_T = BTE.dynamics.material_0.position_0_0_0.electrons.band_0.total.temperature.temperature.nxvalue[BTE_idx]
    BTE_neqdis=BTE_data[BTE_idx,:] - FermiDirac(BTE_energy,BTE_cp,BTE_T)

    AthEM_cp = AthEM.dynamics.material_0.position_0_0_0.electrons.band_0.total.abs_chemical_potential.abs_chemical_potential.nxvalue[AthEM_idx]
    AthEM_T = AthEM.postprocessing.temperature.temperature.nxvalue[AthEM_idx]
    AthEM_ehpdis=AthEM_data[AthEM_idx,:] - FermiDirac(AthEM_energy,AthEM_cp,AthEM_T)

    AthEM_ehp = AthEM_DOS*AthEM_ehpdis
    BTE_ehp = BTE_DOS*BTE_neqdis

    ax.clear()
    ax.text(-4,0.11,str(round(sim_t/1e-15,2))+' fs',fontsize=20)
    ax.plot(BTE_energy[200:-100],BTE_ehp[200:-100],label='Boltzmann',linewidth=5,color='red')
    ax.plot(AthEM_energy[200:-100],AthEM_ehp[200:-100],label='AthEM',linewidth=5,color='blue',dashes=[6,2])
    ax.set_xlabel(r'$\mathrm{E-E_F / eV}$',fontsize=24)
    ax.set_ylabel(r'$\mathrm{Particle \ Density / atom^{-1}}$',fontsize=24)
    ax.set_title("\n".join(wrap('Number of non-equilibrium electron-hole pairs',30)),fontsize=28,fontweight='bold')
    ax.tick_params(axis='both',which='major',labelsize=20,width=3)
    ax.set_ylim(-0.1,0.1)
    ax.legend(prop={'size': 20},loc="upper right")

    return mplfig_to_npimage(fig)

fig,ax=plt.subplots(tight_layout=True)
duration = 12
animation=VideoClip(make_frame,duration=duration)

animation.write_videofile("ehp.mp4",fps=10,codec='libx264',audio=False,threads=8,logger="bar")