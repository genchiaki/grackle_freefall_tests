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
from scipy import interpolate
from yt.units import mp
#from unyt import mp, kboltz, G, h, eV, yr, Myr, Gyr, erg, s, Msun, cm, pc
#from math import pi
from yt.utilities.physical_constants import \
    gravitational_constant_cgs as G

import matplotlib.font_manager as font_manager
font_dirs = ['/mnt/c/Windows/Fonts', ]
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
font_list = font_manager.createFontList(font_files)
font_manager.fontManager.ttflist.extend(font_list)
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams["mathtext.fontset"] = "stix"
fontsize_suptitle = 28
fontsize_title    = 24
fontsize_boxsize  = 24
fontsize_label    = 24
fontsize_cblabel  = 24
fontsize_tick     = 24
fontsize_label_s  = 16
fontsize_cblabel_s= 16
fontsize_tick_s   = 16
fontsize_legend   = 16
fontsize_legend_s = 12

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

i_nH       =  1 - 1
i_Tg       =  2 - 1
i_elec     =  3 - 1
i_HI       =  4 - 1
i_HII      =  5 - 1
i_H2I      =  6 - 1
i_HM       =  7 - 1
i_H2II     =  8 - 1
i_HeHII    =  9 - 1
i_HeI      = 10 - 1
i_HeII     = 11 - 1
i_HeIII    = 12 - 1
i_DI       = 13 - 1
i_DII      = 14 - 1
i_DM       = 15 - 1
i_HDI      = 16 - 1
i_HDII     = 17 - 1
i_CI       = 18 - 1
i_CII      = 19 - 1
i_CO       = 20 - 1
i_CO2      = 21 - 1
i_OI       = 22 - 1
i_OH       = 23 - 1
i_H2O      = 24 - 1
i_O2       = 25 - 1
i_SiI      = 26 - 1
i_SiOI     = 27 - 1
i_SiO2I    = 28 - 1
i_CH       = 29 - 1
i_CH2      = 30 - 1
i_COII     = 31 - 1
i_OII      = 32 - 1
i_OHII     = 33 - 1
i_H2OII    = 34 - 1
i_H3OII    = 35 - 1
i_O2II     = 36 - 1
i_Mg       = 37 - 1
i_Al       = 38 - 1
i_S        = 39 - 1
i_Fe       = 40 - 1
i_SiM      = 41 - 1
i_FeM      = 42 - 1
i_Mg2SiO4  = 43 - 1
i_MgSiO3   = 44 - 1
i_Fe3O4    = 45 - 1
i_AC       = 46 - 1
i_SiO2D    = 47 - 1
i_MgO      = 48 - 1
i_FeS      = 49 - 1
i_Al2O3    = 50 - 1
i_reforg   = 51 - 1
i_volorg   = 52 - 1
i_H2Oice   = 53 - 1
i_tdust    = 54 - 1
i_tSiM     = 55 - 1
i_tFeM     = 56 - 1
i_tMg2SiO4 = 57 - 1
i_tMgSiO3  = 58 - 1
i_tFe3O4   = 69 - 1
i_tAC      = 60 - 1
i_tSiO2D   = 61 - 1
i_tMgO     = 62 - 1
i_tFeS     = 63 - 1
i_tAl2O3   = 64 - 1
i_treforg  = 65 - 1
i_tvolorg  = 66 - 1
i_tH2Oice  = 67 - 1
i_Gadia    = 68 - 1
indir = 'fig'


Zcolors = [ 'red'    , 'orange'    , 'yellow'   , 'lime'     , 'green'
          , 'blue'   , 'purple'    , 'black'    , 'black'    , 'black',   'black']
Zlabels = [
  r'$1$ Z$_{\bigodot}$'      
, r'$10^{-1}$ Z$_{\bigodot}$'
, r'$10^{-2}$ Z$_{\bigodot}$'
, r'$10^{-3}$ Z$_{\bigodot}$'
, r'$10^{-4}$ Z$_{\bigodot}$'
, r'$10^{-5}$ Z$_{\bigodot}$'
, r'$10^{-6}$ Z$_{\bigodot}$'
, r'$10^{-7}$ Z$_{\bigodot}$'
, r'$10^{-8}$ Z$_{\bigodot}$'
, r'$0$ Z$_{\bigodot}$'
         ]
Dcolors = ['red', 'magenta', 'orange', 'green', 'black', 'black', 'black', 'black', 'black', 'blue']
Dlabels = [
  r'$1 G_0$'
, r'$10^{-1} G_0$'
, r'$10^{-2} G_0$'
, r'$10^{-3} G_0$'
, r'$10^{-4} G_0$'
, r'$10^{-5} G_0$'
, r'$10^{-6} G_0$'
, r'$10^{-7} G_0$'
, r'$10^{-8} G_0$'
, r'$0 G_0$'
         ]
SNcolors = ['black'
,'red', 'orange', 'yellow', 'green'
,'black', 'red', 'blue', 'green'
,'blue', 'purple', 'black']
SNlabels = [
  r'Local ISM dust model'
, r'CCSN $M_{\rm PopIII} =  13$ M$_{\bigodot}$'
, r'CCSN $M_{\rm PopIII} =  20$ M$_{\bigodot}$'
, r'CCSN $M_{\rm PopIII} =  25$ M$_{\bigodot}$'
, r'CCSN $M_{\rm PopIII} =  30$ M$_{\bigodot}$'
, r'FSN  $M_{\rm PopIII} =  13$ M$_{\bigodot}$'
, r'FSN  $M_{\rm PopIII} =  15$ M$_{\bigodot}$'
, r'FSN  $M_{\rm PopIII} =  50$ M$_{\bigodot}$'
, r'FSN  $M_{\rm PopIII} =  80$ M$_{\bigodot}$'
, r'PISN $M_{\rm PopIII} = 170$ M$_{\bigodot}$'
, r'PISN $M_{\rm PopIII} = 200$ M$_{\bigodot}$'
, r'Simple dust model'
         ]
GGlinestyles = ['--', '-']

def rho(nH):
    return nH * mH / XH


PLOT_O05        = True
PLOT_O12        = True
PLOT_M14        = True
PLOT_C15        = True
PLOT_P22        = True
PLOT_multimetal = True

data = {}

test = 'O05'
data[test] = np.empty(10, dtype=list)
nZs = [9, 6, 5, 4, 3, 2, 1, 0]
for nZ in nZs:
    fn = test + '/' + 'output_' + 'Z-%d' % (nZ)
    data[test][nZ] = np.loadtxt(fn)

test = 'O12'
data[test] = np.empty([10, 10], dtype=list)
nZs = [9, 5, 4, 3, 2, 1]
nDs = [0, 1, 2, 3, 9]
for nZ in nZs:
    for nD in nDs:
        fn = test + '/' + 'output_' + 'Z-%d' % (nZ) + '_' + 'D-%d' % (nD)
        data[test][nZ, nD] = np.loadtxt(fn)

test = 'M14'
data[test] = np.empty([12, 10, 2], dtype=list)
nSNs = [6, 7, 8]
nZs = [3]
nGGs = [1]
for nSN in nSNs:
  for nZ in nZs:
    for nGG in nGGs:
        fn = test + '/' + 'output_' + 'SN%02d' % (nSN) + '_' + 'Z-%d' % (nZ) + '_' + 'gg%d' % (nGG)
        data[test][nSN, nZ, nGG] = np.loadtxt(fn)

test = 'C15'
data[test] = np.empty([12, 10, 2], dtype=list)
nSNs = [4]
nZs = [6, 5, 4, 3]
nGGs = [0, 1]
for nSN in nSNs:
  for nZ in nZs:
    for nGG in nGGs:
        fn = test + '/' + 'output_' + 'SN%02d' % (nSN) + '_' + 'Z-%d' % (nZ) + '_' + 'gg%d' % (nGG)
        data[test][nSN, nZ, nGG] = np.loadtxt(fn)

test = 'P22'
data[test] = np.empty(10, dtype=list)
nZs = [3]
for nZ in nZs:
    fn = test + '/' + 'output_' + 'Z-%d' % (nZ)
    data[test][nZ] = np.loadtxt(fn)

test = 'multimetal'
data[test] = np.empty(10, dtype=list)
nZs = [3]
for nZ in nZs:
    fn = test + '/' + 'output_' + 'Z-%d' % (nZ)
    data[test][nZ] = np.loadtxt(fn)


if PLOT_O05:
    test = 'O05'
    nZs = [9, 6, 5, 4, 3, 2, 1, 0]

    nrows = 1
    ncols = 1
    mergin_bottom = 0.9
    mergin_top    = 0.6
    mergin_left   = 0.8
    mergin_right  = 0.2
    figsize_x = mergin_left + ncols*6 + mergin_right
    figsize_y = mergin_bottom + nrows*4 + mergin_top
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(figsize_x, figsize_y))
    fig.suptitle('Omukai et al. (2005)', fontsize=fontsize_suptitle)
    fig.subplots_adjust(left=mergin_left/figsize_x, right=1-mergin_right/figsize_x
                      , bottom=mergin_bottom/figsize_y, top=1-mergin_top/figsize_y
                      , wspace=0.0, hspace=0.0)
    for irow in range(nrows):
      for icol in range(ncols):
        ax = axs
        for nZ in nZs:
            ax.plot(np.log10(data[test][nZ][:,i_nH]), np.log10(data[test][nZ][:,i_Tg])
                 , color=Zcolors[nZ], linestyle='-' , linewidth=2, label=Zlabels[nZ])
            ax.plot(np.log10(data[test][nZ][:,i_nH]), np.log10(data[test][nZ][:,i_tdust])
                 , color=Zcolors[nZ], linestyle='--', linewidth=2, label=None)
        ax.set_xticks(np.linspace( 0, 16, 5))
        ax.set_xticks(np.linspace(-1, 16, 18), minor=True)
        ax.set_yticks(np.linspace( 1, 4, 4))
        ax.set_yticks(np.linspace( 0, 4, 9), minor=True)
        ax.set_xlim(-1, 16)
        ax.set_ylim(0, 4)
        ax.tick_params(labelsize=fontsize_tick)
        if irow < nrows - 1:
            ax.set_xlabel('')
            ax.tick_params(labelbottom = False)
        else:
            ax.set_xlabel(r'log [ Density / cm$^{-3}$ ]', fontsize=fontsize_label)
            ax.tick_params(labelbottom = True)
        if icol > 0:
            ax.set_ylabel('')
            ax.tick_params(labelleft = False)
        else:
            ax.set_ylabel(r'log [ Temperature / K ]'    , fontsize=fontsize_tick)
            ax.tick_params(labelleft = True)
        if irow == nrows - 1 and icol == 0:
            ax.legend(fontsize=fontsize_legend_s, labelspacing=0.1, ncol=2)
    fig.savefig("%s/%s_nT.png" % (indir, test))
    plt.close('all')


if PLOT_O12:
    test = 'O12'
    nZs = [9, 5, 4, 3, 2, 1]
    nDs = [0, 1, 2, 3, 9]
   
    nrows = 3
    ncols = 2
    mergin_bottom = 0.9
    mergin_top    = 0.6
    mergin_left   = 0.8
    mergin_right  = 0.2
    figsize_x = mergin_left + ncols*5 + mergin_right
    figsize_y = mergin_bottom + nrows*4 + mergin_top
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(figsize_x, figsize_y))
    fig.suptitle('Omukai et al. (2012)', fontsize=fontsize_suptitle)
    fig.subplots_adjust(left=mergin_left/figsize_x, right=1-mergin_right/figsize_x
                      , bottom=mergin_bottom/figsize_y, top=1-mergin_top/figsize_y
                      , wspace=0.0, hspace=0.0)
    for irow in range(nrows):
      for icol in range(ncols):
        ax = axs[irow][icol]
        nZ = nZs[icol + ncols * irow]
        ax.annotate(Zlabels[nZ], xy=(0.5, 0.98), xycoords=ax.transAxes
              , color='black', size=fontsize_title, rotation=0
              , va='top', ha='center'
                )
        for nD in nDs:
            ax.plot(np.log10(data[test][nZ,nD][:,i_nH]), np.log10(data[test][nZ,nD][:,i_Tg]), color=Dcolors[nD], linestyle='-' , linewidth=2, label=Dlabels[nD])
        ax.set_xticks(np.linspace( 0, 16, 5))
        ax.set_xticks(np.linspace(-1, 16, 18), minor=True)
        ax.set_yticks(np.linspace( 1, 4, 4))
        ax.set_yticks(np.linspace( 0, 4, 9), minor=True)
        ax.set_xlim(-1, 12)
        ax.set_ylim(1, 4.5)
        ax.tick_params(labelsize=fontsize_tick)
        if irow < nrows - 1:
            ax.set_xlabel('')
            ax.tick_params(labelbottom = False)
        else:
            ax.set_xlabel(r'log [ Density / cm$^{-3}$ ]', fontsize=fontsize_label)
            ax.tick_params(labelbottom = True)
        if icol > 0:
            ax.set_ylabel('')
            ax.tick_params(labelleft = False)
        else:
            ax.set_ylabel(r'log [ Temperature / K ]'    , fontsize=fontsize_tick)
            ax.tick_params(labelleft = True)
        if irow == nrows - 1 and icol == 0:
            ax.legend(fontsize=fontsize_legend, labelspacing=0.1, ncol=1)
    fig.savefig("%s/%s_nT.png" % (indir, test))
    plt.close('all')


if PLOT_M14:
    test = 'M14'
    nZs = [3]
    nSNs = [6, 7, 8]
    nGGs = [1]
   
    nrows = 1
    ncols = 3
    mergin_bottom = 0.9
    mergin_top    = 0.6
    mergin_left   = 0.8
    mergin_right  = 0.2
    figsize_x = mergin_left + ncols*5 + mergin_right
    figsize_y = mergin_bottom + nrows*4 + mergin_top
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(figsize_x, figsize_y))
    fig.suptitle('Marassi et al. (2014)', fontsize=fontsize_suptitle)
    fig.subplots_adjust(left=mergin_left/figsize_x, right=1-mergin_right/figsize_x
                      , bottom=mergin_bottom/figsize_y, top=1-mergin_top/figsize_y
                      , wspace=0.0, hspace=0.0)
    for irow in range(nrows):
      for icol in range(ncols):
        ax = axs[icol]
        nSN = nSNs[icol + ncols * irow]
        nZ  = 3
        nGG = 1
        ax.annotate(SNlabels[nSN], xy=(0.5, 0.98), xycoords=ax.transAxes
              , color='black', size=fontsize_title, rotation=0
              , va='top', ha='center'
                )
        ax.plot(np.log10(data[test][nSN,nZ,nGG][:,i_nH]), np.log10(data[test][nSN,nZ,nGG][:,i_Tg])
             , color=SNcolors[nSN], linestyle='-' , linewidth=2, label=SNlabels[nSN])
        ax.set_xticks(np.linspace( 0, 16, 5))
        ax.set_xticks(np.linspace(-1, 16, 18), minor=True)
        ax.set_yticks(np.linspace( 1, 4, 4))
        ax.set_yticks(np.linspace( 0, 4, 9), minor=True)
        ax.set_xlim(-1, 16)
        ax.set_ylim(1, 4)
        ax.tick_params(labelsize=fontsize_tick)
        if irow < nrows - 1:
            ax.set_xlabel('')
            ax.tick_params(labelbottom = False)
        else:
            ax.set_xlabel(r'log [ Density / cm$^{-3}$ ]', fontsize=fontsize_label)
            ax.tick_params(labelbottom = True)
        if icol > 0:
            ax.set_ylabel('')
            ax.tick_params(labelleft = False)
        else:
            ax.set_ylabel(r'log [ Temperature / K ]'    , fontsize=fontsize_tick)
            ax.tick_params(labelleft = True)
    fig.savefig("%s/%s_nT.png" % (indir, test))
    plt.close('all')


if PLOT_C15:
    test = 'C15'
    nSNs = [4]
    nZs = [6, 5, 4]
    nGGs = [0, 1]

    nrows = 1
    ncols = 1
    mergin_bottom = 0.9
    mergin_top    = 0.6
    mergin_left   = 0.8
    mergin_right  = 0.2
    figsize_x = mergin_left + ncols*5 + mergin_right
    figsize_y = mergin_bottom + nrows*4 + mergin_top
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(figsize_x, figsize_y))
    fig.suptitle('Chiaki et al. (2015)', fontsize=fontsize_suptitle)
    fig.subplots_adjust(left=mergin_left/figsize_x, right=1-mergin_right/figsize_x
                      , bottom=mergin_bottom/figsize_y, top=1-mergin_top/figsize_y
                      , wspace=0.0, hspace=0.0)
    for irow in range(nrows):
      for icol in range(ncols):
        ax = axs
        nSN = 4
        ax.annotate(SNlabels[nSN], xy=(0.5, 0.98), xycoords=ax.transAxes
              , color='black', size=fontsize_title, rotation=0
              , va='top', ha='center'
                )
        for nZ in nZs:
            for nGG in nGGs:
                ax.plot(np.log10(data[test][nSN,nZ,nGG][:,i_nH]), np.log10(data[test][nSN,nZ,nGG][:,i_Tg])
                     , color=Zcolors[nZ], linestyle=GGlinestyles[nGG], linewidth=2, label=(Zlabels[nZ] if nGG else None))
        ax.set_xticks(np.linspace( 0, 16, 5))
        ax.set_xticks(np.linspace(-1, 16, 18), minor=True)
        ax.set_yticks(np.linspace( 1, 4, 4))
        ax.set_yticks(np.linspace( 0, 4, 9), minor=True)
        ax.set_xlim(-1, 16)
        ax.set_ylim(1.5, 3.5)
        ax.tick_params(labelsize=fontsize_tick)
        if irow < nrows - 1:
            ax.set_xlabel('')
            ax.tick_params(labelbottom = False)
        else:
            ax.set_xlabel(r'log [ Density / cm$^{-3}$ ]', fontsize=fontsize_label)
            ax.tick_params(labelbottom = True)
        if icol > 0:
            ax.set_ylabel('')
            ax.tick_params(labelleft = False)
        else:
            ax.set_ylabel(r'log [ Temperature / K ]'    , fontsize=fontsize_tick)
            ax.tick_params(labelleft = True)
        if irow == nrows - 1 and icol == 0:
            ax.legend(fontsize=fontsize_legend, labelspacing=0.1, ncol=1)
    fig.savefig("%s/%s_nT.png" % (indir, test))
    plt.close('all')


if PLOT_P22:
    test = 'P22'
    nZs = [3]

    nrows = 1
    ncols = 1
    mergin_bottom = 0.9
    mergin_top    = 0.6
    mergin_left   = 0.8
    mergin_right  = 0.2
    mergin_h = 0.1
    mergin_v = 0.0
    figsize_x = mergin_left + ncols*5 + (ncols-1)*mergin_h + mergin_right
    figsize_y = mergin_bottom + nrows*4 + (nrows-1)*mergin_v + mergin_top
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(figsize_x, figsize_y))
    fig.suptitle('Park et al. (2022)', fontsize=fontsize_suptitle)
    fig.subplots_adjust(left=mergin_left/figsize_x, right=1-mergin_right/figsize_x
                      , bottom=mergin_bottom/figsize_y, top=1-mergin_top/figsize_y
                      , wspace=mergin_h, hspace=mergin_v)
    for irow in range(nrows):
      for icol in range(ncols):
        ax = axs # [icol]
        nZ = 3
        if irow == 0 and icol == 0:
            ax.plot(np.log10(data[test][nZ][:,i_nH]), np.log10(data[test][nZ][:,i_Tg])
                 , color='blue' , linestyle='-' , linewidth=2, label=r'$T_{\rm gas}$')
##          ax.plot(np.log10(data[test][nZ][:,i_nH]), np.log10(data[test][nZ][:,i_tMgSiO3])
##               , color='green', linestyle='--', linewidth=2, label=r'$T_{\rm silicate}$')
            ax.plot(np.log10(data[test][nZ][:,i_nH]), np.log10(data[test][nZ][:,i_tAC])
                 , color='black', linestyle='--', linewidth=2, label=r'$T_{\rm graphite}$')
        if irow == 0 and icol == 1:
            ax.plot(np.log10(data[test][nZ][:,i_nH]), np.log10(data[test][nZ][:,i_MgSiO3])
                 , color='green', linestyle='-', linewidth=2, label=r'$silicate$')
            ax.plot(np.log10(data[test][nZ][:,i_nH]), np.log10(data[test][nZ][:,i_AC])
                 , color='black', linestyle='-', linewidth=2, label=r'$graphite$')
        ax.set_xticks(np.linspace( 0, 16, 5))
        ax.set_xticks(np.linspace(-1, 16, 18), minor=True)
        ax.set_xlim(-1, 16)
        if irow == 0 and icol == 0:
            ax.set_yticks(np.linspace( 1, 4, 4))
            ax.set_yticks(np.linspace( 0, 4, 9), minor=True)
            ax.set_ylim(1.5, 3.5)
        ax.tick_params(labelsize=fontsize_tick)
        if irow < nrows - 1:
            ax.set_xlabel('')
            ax.tick_params(labelbottom = False)
        else:
            ax.set_xlabel(r'log [ Density / cm$^{-3}$ ]', fontsize=fontsize_label)
            ax.tick_params(labelbottom = True)
        if irow == 0 and icol == 0:
            ax.set_ylabel(r'log [ Temperature / K ]'    , fontsize=fontsize_tick)
        if irow == 0 and icol == 1:
            ax.set_ylabel(r'log [ Abundance ]'         , fontsize=fontsize_tick)
        ax.tick_params(labelleft = True)
        if irow == nrows - 1 and icol == 0:
            ax.legend(fontsize=fontsize_legend, labelspacing=0.1, ncol=1)
    fig.savefig("%s/%s_nT.png" % (indir, test))
    plt.close('all')


if PLOT_multimetal:
    test = 'multimetal'
    nZs = [3]

    nrows = 1
    ncols = 1
    mergin_bottom = 0.9
    mergin_top    = 0.6
    mergin_left   = 0.8
    mergin_right  = 0.2
    figsize_x = mergin_left + ncols*5 + mergin_right
    figsize_y = mergin_bottom + nrows*4 + mergin_top
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(figsize_x, figsize_y))
    fig.suptitle('MultiMetals', fontsize=fontsize_suptitle)
    fig.subplots_adjust(left=mergin_left/figsize_x, right=1-mergin_right/figsize_x
                      , bottom=mergin_bottom/figsize_y, top=1-mergin_top/figsize_y
                      , wspace=0.0, hspace=0.0)
    for irow in range(nrows):
      for icol in range(ncols):
        ax = axs
        test = 'C15'; nSN = 4; nZ = 3; nGG = 1
        ax.plot(np.log10(data[test][nSN,nZ,nGG][:,i_nH]), np.log10(data[test][nSN,nZ,nGG][:,i_Tg])
             , color='blue', linestyle='--', linewidth=2, label=SNlabels[nSN]+' (C30)')
        test = 'M14'; nSN = 6; nZ = 3; nGG = 1
        ax.plot(np.log10(data[test][nSN,nZ,nGG][:,i_nH]), np.log10(data[test][nSN,nZ,nGG][:,i_Tg])
             , color='red', linestyle='--', linewidth=2, label=SNlabels[nSN]+' (F15)')
        test = 'multimetal'; nZ = 3
        ax.plot(np.log10(data[test][nZ][:,i_nH]), np.log10(data[test][nZ][:,i_Tg])
             , color='purple', linestyle='-', linewidth=2, label='50% C30 + 50% F15')
        ax.set_xticks(np.linspace( 0, 16, 5))
        ax.set_xticks(np.linspace(-1, 16, 18), minor=True)
        ax.set_yticks(np.linspace( 1, 4, 4))
        ax.set_yticks(np.linspace( 0, 4, 9), minor=True)
        ax.set_xlim(-1, 16)
        ax.set_ylim(1.5, 3.5)
        ax.tick_params(labelsize=fontsize_tick)
        if irow < nrows - 1:
            ax.set_xlabel('')
            ax.tick_params(labelbottom = False)
        else:
            ax.set_xlabel(r'log [ Density / cm$^{-3}$ ]', fontsize=fontsize_label)
            ax.tick_params(labelbottom = True)
        if icol > 0:
            ax.set_ylabel('')
            ax.tick_params(labelleft = False)
        else:
            ax.set_ylabel(r'log [ Temperature / K ]'    , fontsize=fontsize_tick)
            ax.tick_params(labelleft = True)
        if irow == nrows - 1 and icol == 0:
            ax.legend(fontsize=fontsize_legend, labelspacing=0.1, ncol=1)
    fig.savefig("%s/%s_nT.png" % (indir, test))
    plt.close('all')


