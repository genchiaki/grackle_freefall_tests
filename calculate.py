import yt
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
#from yt.data_objects.particle_filters import add_particle_filter
#from yt.analysis_modules.halo_finding.api import HaloFinder
from yt.extensions.astro_analysis.halo_analysis.halo_catalog import HaloCatalog
from yt.data_objects.particle_filters import add_particle_filter
from yt.extensions.astro_analysis.halo_finding.rockstar.api import RockstarHaloFinder
'''
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
from yt.data_objects.particle_filters import add_particle_filter
from yt.extensions.astro_analysis.halo_analysis.halo_finding.rockstar.api import RockstarHaloFinder
'''
import numpy as np
from numpy import linalg as LA
import datetime
import sys
import struct
from yt.units import mp
#from unyt import mp, kboltz, G, h, eV, yr, Myr, Gyr, erg, s, Msun, cm, pc
#from math import pi
from yt.utilities.physical_constants import \
    gravitational_constant_cgs as G

pi = 3.14159265
H0 = 100.0 * 1e5 / 3.0856e24 # cgs
Delta = 18.0 * pi**2
boltzman = 1.38e-16
Msun = 1.989e33
HubbleParam = 0.71
Omega0      = 0.314085 
OmegaLambda = 0.685915 
f_baryon    = 0.0440119 / Omega0
GRAVITY = 6.672e-8
mH = 1.6726e-24
mass_range = (7.0,10.0)
boxsize = 0.6774 # Mpc
max_spin = 0.06
dtDataDump0 = 0.81650437825362
dtDataDump = 0.733626
TimeUnits = 3.996888808944e+14
Zsun = 0.02
yr_s  = 3.1556952e7
Myr_s = 3.1556952e13
XH = 0.76
au_cm = 1.49598e13
pc_cm = 3.0856e18
sigma_sb = 5.670373e-5

HydrogenFractionByMass   = 0.76
DeuteriumToHydrogenRatio = 3.4e-5 * 2.0
HeliumToHydrogenRatio    = (1.0 - HydrogenFractionByMass) / HydrogenFractionByMass
SolarMetalFractionByMass = 0.01295
SolarIronAbundance = 7.50

UnitLength_in_cm= 3.085678e21
UnitMass_in_g= 1.989e43
UnitVelocity_in_cm_per_s= 1.0e5
UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s
UnitDensity_in_cgs= UnitMass_in_g/ (UnitLength_in_cm**3)
UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ (UnitTime_in_s**2)
UnitEnergy_in_cgs= UnitMass_in_g * (UnitLength_in_cm**2) / (UnitTime_in_s**2)
G=GRAVITY/ (UnitLength_in_cm**3) * UnitMass_in_g * (UnitTime_in_s**2)

MinimumComovingHydroSoftening = 0.03
GasSoftFactor                 = 2.5
ReferenceGasPartMass          = 3.9012e-10

DUST_TEMPERATURE = True

def rho2nH(time, rho):
    return (rho * HubbleParam**2 / time**3 * UnitDensity_in_cgs) * HydrogenFractionByMass / mH

def rho(nH):
    return nH * mH / XH

if DUST_TEMPERATURE:
  
    # isrf ################################################
    # in cool1d_multi_g.F
    fgr = 0.009387 # dust to gas mass ratio
    isrf = 2.00349e12 # 1e3 # habing
    sigma_MgSiO3 = 7.5835e-4 # cm^2
    gamma_isrf_MgSiO3 = 5.3e-3 * sigma_MgSiO3
    # in calc_tdust_1d_g.F
    gamma_isrf_MgSiO3 *= isrf
    print(gamma_isrf_MgSiO3)

    # radiation from grains ###############################
    Trad = (55.545718409999992)**0.25
    TMgSiO3 = 299.72521132488612 # 52.947188743978494
    kappa_MgSiO3 = 2.8612213552560130E-005 # 2.2547882054727010E-006 # cm^2 / gas g
    gamma_rad_MgSiO3 = 4.0 * sigma_sb * kappa_MgSiO3 * (Trad**4 - TMgSiO3**4)
    print(4.0 * sigma_sb, kappa_MgSiO3)
    print(Trad**4, TMgSiO3)
#   print(gamma_rad_MgSiO3)

    # heat transfer from gas to dust ######################
    # in cool1d_multi_g.F
    units = 3.21310e-28
    nh = 0.1
    Tgas = 300
    f_vel = 0.5 / np.sqrt(2.0) + 0.0833333 / np.sqrt(4.0)
    units = 1.0
    vH_avg = np.sqrt(boltzman * Tgas / 2.0 / pi / mH)
    fv2k = f_vel * 4.0 * vH_avg * 2.0 * boltzman # * mH / units
#   fac = units / mH
    gasMgSiO3 = fv2k * sigma_MgSiO3 # * fac
    print(gasMgSiO3)
    # in calc_tdust_1d_g.F
    gamma_gas_MgSiO3 = gasMgSiO3 * nh * (Tgas - TMgSiO3)
    print(gamma_gas_MgSiO3)

#   print(gamma_isrf_MgSiO3 + gamma_rad_MgSiO3 + gamma_gas_MgSiO3)
