#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""HMF plots"""

import functools
import numpy as np
import os

import common
import utilities_statistics as us

##################################

# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3
c_light = 299792458.0 #m/s

# Mass function initialization

mlow = -30 + 5.0 * np.log10(0.677)
mupp = -10 + 5.0 * np.log10(0.677)
dm = 0.5
mbins = np.arange(mlow,mupp,dm)
xlf   = mbins + dm/2.0


mflow = 5
mfupp = 14
dmf = 0.5
mfbins = np.arange(mflow,mfupp,dmf)
xmf = mfbins + dmf/2.0

mflow2 = 9
mfupp2 = 14
dmf2 = 0.25
mfbins2 = np.arange(mflow2,mfupp2,dmf2)
xmf2 = mfbins2 + dmf2/2.0


def plot_uv_mass_microsurfs(x, z, color = 'black', stellar_mass = True, scatter=True, label = True, labelin="z=6"):

    if(stellar_mass == True):
       sm_micro,luv_micro,suvdn_micro,suvup_micro,z_micro = np.loadtxt("../data/Models/SharkVariations/uv_stellarmass_L24HBTbestparams_micro-SURFS.txt", unpack = True, usecols = [0,1,2,3,4])
    else:
       sm_micro,luv_micro,suvdn_micro,suvup_micro,z_micro = np.loadtxt("../data/Models/SharkVariations/uv_halomass_L24HBTbestparams_micro-SURFS.txt", unpack = True, usecols = [0,1,2,3,4])
    
    ind=np.where((z_micro == z) & (luv_micro != 0))
    if(scatter == False):
        ax.plot(sm_micro[ind], luv_micro[ind], color=color,linestyle='solid', linewidth = 1) #, label = labelin if label == True else 'None')
    else:
        ax.plot(sm_micro[ind], suvdn_micro[ind], color=color, linestyle='solid', linewidth = 1)
        ax.plot(sm_micro[ind], suvup_micro[ind], color=color, linestyle='solid', linewidth = 1)


def plot_uv_mass_l800(x, z, color = 'black', stellar_mass = True, scatter=True, label = True, labelin="z=6", uvlim = -13):

    if(stellar_mass == True):
       sm_l800,luv_l800,suvdn_l800,suvup_l800,z_l800 = np.loadtxt("../data/Models/SharkVariations/uv_stellarmass_inverse_L24HBTbestparams_L800.txt", unpack = True, usecols = [0,1,2,3,4])
    else:
       sm_l800,luv_l800,suvdn_l800,suvup_l800,z_l800 = np.loadtxt("../data/Models/SharkVariations/uv_halomass_inverse_L24HBTbestparams_L800.txt", unpack = True, usecols = [0,1,2,3,4])

    ind=np.where((z_l800 == z) & (luv_l800 != 0) & (sm_l800 <= uvlim))
    if(scatter == False):
        ax.plot(sm_l800[ind], luv_l800[ind], color=color,linestyle='solid', linewidth = 1.5) #, label = labelin if label == True else 'None')
    else:
        ax.plot(sm_l800[ind], -1 * suvdn_l800[ind], color=color, linestyle='solid', linewidth = 1.5)
        ax.plot(sm_l800[ind], suvup_l800[ind], color=color, linestyle='solid', linewidth = 1.5)

def plot_uv_mass_l800_nodiskins(x, z, color = 'black', stellar_mass = True, scatter=True, label = True, labelin="z=6", uvlim = -13):

    if(stellar_mass == True):
       sm_l800,luv_l800,suvdn_l800,suvup_l800,z_l800 = np.loadtxt("../data/Models/SharkVariations/uv_stellarmass_inverse_L24HBTbestparams_nodiskins_L800.txt", unpack = True, usecols = [0,1,2,3,4])
    else:
       sm_l800,luv_l800,suvdn_l800,suvup_l800,z_l800 = np.loadtxt("../data/Models/SharkVariations/uv_halomass_inverse_L24HBTbestparams_nodiskins_L800.txt", unpack = True, usecols = [0,1,2,3,4])
    #suvd_l800 = (abs(suvdn_l800) + abs(suvup_l800) ) /2.   
    ind=np.where((z_l800 == z) & (luv_l800 != 0) & (sm_l800 <= uvlim))
    if(scatter == False):
        ax.plot(sm_l800[ind], luv_l800[ind], color=color,linestyle='dashed', linewidth = 1) #, label = labelin if label == True else 'None')
    else:
        ax.plot(sm_l800[ind], -1. * suvdn_l800[ind], color=color, linestyle='dashed', linewidth = 1)
        ax.plot(sm_l800[ind], suvup_l800[ind], color=color, linestyle='dashed', linewidth = 1)


def plot_uv_mass_hisurfs(x, z, color = 'black', stellar_mass = True, scatter=True, label = True, labelin="z=6"):

    if(stellar_mass == True):
       sm_micro,luv_micro,suvdn_micro,suvup_micro,z_micro = np.loadtxt("../data/Models/SharkVariations/uv_stellarmass_L24HBTbestparams_HiSURFS_L50.txt", unpack = True, usecols = [0,1,2,3,4])
    else:
       sm_micro,luv_micro,suvdn_micro,suvup_micro,z_micro = np.loadtxt("../data/Models/SharkVariations/uv_halomass_L24HBTbestparams_HiSURFS_L50.txt", unpack = True, usecols = [0,1,2,3,4])
    
    ind=np.where((z_micro == z) & (luv_micro != 0))
    if(scatter == False):
        ax.plot(sm_micro[ind], luv_micro[ind], color=color,linestyle='dotted', linewidth = 1) #, label = labelin if label == True else 'None')
    else:
        ax.plot(sm_micro[ind], -1*suvdn_micro[ind], color=color, linestyle='dotted', linewidth = 1)
        ax.plot(sm_micro[ind], suvup_micro[ind], color=color, linestyle='dotted', linewidth = 1)


def plot_uv_contributions_l800(ax, z, label = True):
    x, ysfd, ydi, ym, redshift = np.loadtxt("../data/Models/SharkVariations/uv_mag_channels_L24HBTbestparams_L800_allz.txt", unpack = True, usecols = [0,1,2,3,4])

    ind = np.where((redshift == z) & (ysfd != -1))
    ylow = np.zeros(shape = len(x[ind]))
    ax.fill_between(x[ind], ylow, ysfd[ind], facecolor='navy', alpha=0.15, interpolate=True, label = 'SF disks (L800)' if label == True else None)
    ax.fill_between(x[ind], ysfd[ind], ysfd[ind] + ydi[ind], facecolor='DarkSeaGreen', alpha=0.15, interpolate=True, label = 'SBs disk ins. (L800)' if label == True else None)
    yhigh = np.zeros(shape = len(x[ind]))
    yhigh[:] = 1.0
    ax.fill_between(x[ind], ysfd[ind] + ydi[ind], yhigh, facecolor='Firebrick', alpha=0.15, interpolate=True, label = 'SBs mergers (L800)' if label == True else None)


def plot_uv_contributions_hisurfs(ax, z, label = True):
    x, ysfd, ydi, ym, redshift = np.loadtxt("../data/Models/SharkVariations/uv_mag_channels_L24HBTbestparams_HiSURFS_L50_allz.txt", unpack = True, usecols = [0,1,2,3,4])

    ind = np.where((redshift == z) & (ysfd != -1))
    ax.plot(x[ind], ysfd[ind], linestyle='solid', color='navy', label = 'SF disks (Hi-SURFS)' if label == True else 'None')
    #ax.plot(x[ind], ysfd[ind] + ydi[ind], linestyle='solid', color='DarkSeaGreen') #, label='SBs disk ins.' if label == True else None)


def plot_uv_contributions_microsurfs(ax, z, label = True):

    x, ysfd, ydi, ym, redshift = np.loadtxt("../data/Models/SharkVariations/uv_mag_channels_L24HBTbestparams_micro-SURFS_allz.txt", unpack = True, usecols = [0,1,2,3,4])

    ind = np.where((redshift == z) & (ysfd != -1))
    ax.plot(x[ind], ysfd[ind], linestyle='dashed', color='navy', label = 'SF disks (micro-SURFS)' if label == True else 'None')
    #ax.plot(x[ind], ysfd[ind] + ydi[ind], linestyle='solid', color='DarkSeaGreen') #, label='SBs disk ins.' if label == True else None)


def plot_uv_lf_hisurfs(ax, z, all_channels = False, label = True):
    x, y, yun, ydi, ym, yd, redshift = np.loadtxt("../data/Models/SharkVariations/uv_lf_L24HBTbestparams_HiSURFS_L50_allz.txt", unpack = True, usecols = [0,1,2,3,4,5,6])

    if(all_channels == False):
       if(z != 17):
          ind = np.where((y < -0.3) & (redshift == z))
          ax.plot(x[ind], y[ind], color='Orchid', linewidth=2, label = 'Total att. (Hi)' if label == True else None)

          ind = np.where((yun < -0.3) & (redshift == z))
          ax.plot(x[ind], yun[ind],color='Orchid', linewidth=2, linestyle='dashed', label='Total unatt. (Hi)' if label == True else None)
       else:
          ind = np.where((y < -0.3) & (redshift == z))
          ax.plot(x[ind], y[ind], color='Orchid', linewidth=1, label = 'Total att. (Hi)' if label == True else None)

    if(all_channels):
       ind = np.where((ydi < -0.3) & (redshift == z))
       ax.plot(x[ind], ydi[ind],'Navy', linewidth=2, linestyle='dashed', label = 'SF disks (Hi)' if label == True else None)
       ind = np.where((ym < -0.3) & (redshift == z))
       ax.plot(x[ind], ym[ind], 'DarkSeaGreen', linewidth=2, linestyle='dashed', label='SBs disk ins. (Hi)' if label == True else None)
       ind = np.where((yd < -0.3) & (redshift == z))
       ax.plot(x[ind], yd[ind],'Firebrick', linewidth=2, linestyle='dashed', label='SBs mergers (Hi)' if label == True else None)
 

def plot_uv_lf_medi(ax, z, all_channels = False, label = True):
    x, y, yun, ydi, ym, yd, redshift = np.loadtxt("../data/Models/SharkVariations/uv_lf_L24HBTbestparams_medi-SURFS_allz.txt", unpack = True, usecols = [0,1,2,3,4,5,6])

    if(all_channels == False):
       if(z != 17):
          ind = np.where((y < -0.3) & (redshift == z))
          ax.plot(x[ind], y[ind], color='black', linewidth=2, label = 'Total att. (medi)' if label == True else None)

          ind = np.where((yun < -0.3) & (redshift == z))
          ax.plot(x[ind], yun[ind],color='black', linewidth=2, linestyle='dashed', label='Total unatt. (medi)' if label == True else None)
       else:
          ind = np.where((y < -0.3) & (redshift == z))
          ax.plot(x[ind], y[ind], color='black', linewidth=1, label = 'Total att. (medi)' if label == True else None)

    if(all_channels):
       ind = np.where((ydi < -0.3) & (redshift == z))
       ax.plot(x[ind], ydi[ind],'Navy', linewidth=2, linestyle='dashed', label = 'SF disks (micro)' if label == True else None)
       ind = np.where((ym < -0.3) & (redshift == z))
       ax.plot(x[ind], ym[ind], 'DarkSeaGreen', linewidth=2, linestyle='dashed', label='SBs disk ins. (micro)' if label == True else None)
       ind = np.where((yd < -0.3) & (redshift == z))
       ax.plot(x[ind], yd[ind],'Firebrick', linewidth=2, linestyle='dashed', label='SBs mergers (micro)' if label == True else None)
  
def plot_uv_lf_micro(ax, z, all_channels = False, label = True):
    x, y, yun, ydi, ym, yd, redshift = np.loadtxt("../data/Models/SharkVariations/uv_lf_L24HBTbestparams_micro-SURFS_allz.txt", unpack = True, usecols = [0,1,2,3,4,5,6])

    if(all_channels == False):
       if(z != 17):
          ind = np.where((y < -0.3) & (redshift == z))
          ax.plot(x[ind], y[ind], color='Purple', linewidth=2, label = 'Total att. (micro)' if label == True else None)

          ind = np.where((yun < -0.3) & (redshift == z))
          ax.plot(x[ind], yun[ind],color='Purple', linewidth=2, linestyle='dashed', label='Total unatt. (micro)' if label == True else None)
       else:
          ind = np.where((y < -0.3) & (redshift == z))
          ax.plot(x[ind], y[ind], color='Purple', linewidth=1, label = 'Total att. (micro)' if label == True else None)

    if(all_channels):
       ind = np.where((ydi < -0.3) & (redshift == z))
       ax.plot(x[ind], ydi[ind],'Navy', linewidth=2, linestyle='dashed', label = 'SF disks (micro)' if label == True else None)
       ind = np.where((ym < -0.3) & (redshift == z))
       ax.plot(x[ind], ym[ind], 'DarkSeaGreen', linewidth=2, linestyle='dashed', label='SBs disk ins. (micro)' if label == True else None)
       ind = np.where((yd < -0.3) & (redshift == z))
       ax.plot(x[ind], yd[ind],'Firebrick', linewidth=2, linestyle='dashed', label='SBs mergers (micro)' if label == True else None)
   
def plot_uv_lf_l800(ax, z, all_channels = False, label = True):
    x, y, yun, ydi, ym, yd, redshift = np.loadtxt("../data/Models/SharkVariations/uv_lf_L24HBTbestparams_L800_allz.txt", unpack = True, usecols = [0,1,2,3,4,5,6])

    if(all_channels == False):
       if(z != 17):
          ind = np.where((y < -0.3) & (redshift == z))
          ax.plot(x[ind], y[ind], color='Olive', linewidth=2, label = 'Total att. (L800)' if label == True else None)

          ind = np.where((yun < -0.3) & (redshift == z))
          ax.plot(x[ind], yun[ind],color='Olive', linewidth=2, linestyle='dashed', label='Total unatt. (L800)' if label == True else None)
       else:
          ind = np.where((y < -0.3) & (redshift == z))
          ax.plot(x[ind], y[ind], color='Olive', linewidth=1, label = 'Total att. (L800)' if label == True else None)


    if(all_channels):
       ind = np.where(ydi != -10)
       ax.plot(x[ind], ydi[ind],'Navy', linewidth=2, linestyle='dotted', label = 'SF disks (L800)' if label == True else None)
       ind = np.where(ym != -10)
       ax.plot(x[ind], ym[ind], 'DarkSeaGreen', linewidth=2, linestyle='dotted', label='SBs disk ins. (L800)' if label == True else None)
       ind = np.where(yd != -10)
       ax.plot(x[ind], yd[ind],'Firebrick', linewidth=2, linestyle='dotted', label='SBs mergers (L800)' if label == True else None)
   

def plot_uv_lf_l800_variants(ax, zl, z, all_channels = False, label = True, color='black'):
    x, y, yun, ydi, ym, yd, zin = np.loadtxt("../data/Models/SharkVariations/uv_lf_L24HBTbestparams_L800_allz.txt", unpack = True, usecols = [0,1,2,3,4,5,6])

    if(all_channels == False):
        ind = np.where((y < -1) & (zin == z))
        ax.plot(x[ind], y[ind], color=c, linewidth=3, label = 'fiducial' if label == True else None)

    x, y, yun, ydi, ym, yd, zin = np.loadtxt("../data/Models/SharkVariations/uv_lf_L24HBTbestparams_nodiskins_L800_allz.txt", unpack = True, usecols = [0,1,2,3,4,5,6])

    if(all_channels == False):
        ind = np.where((y < -1) & (zin == z))
        ax.plot(x[ind], y[ind], color=c, linewidth=3, linestyle='dashed', label = 'no disk ins.' if label == True else None)

    x, y, yun, ydi, ym, yd, zin = np.loadtxt("../data/Models/SharkVariations/uv_lf_L24HBTbestparams_boost1_L800_allz.txt", unpack = True, usecols = [0,1,2,3,4,5,6])

    if(all_channels == False):
        ind = np.where((y < -1) & (zin == z))
        ax.plot(x[ind], y[ind], color=c, linewidth=3, linestyle='dotted', label = '$\\eta_{\\rm burst}=1$' if label == True else None)


def plot_uv_lf_l800_dust_variants(ax, zl, z, all_channels = False, label = True, color='black'):
    x, y, yun, ydi, ym, yd, zin = np.loadtxt("../data/Models/SharkVariations/uv_lf_L24HBTbestparams_L800_allz.txt", unpack = True, usecols = [0,1,2,3,4,5,6])

    if(all_channels == False):
        ind = np.where((y < -1) & (zin == z))
        ax.plot(x[ind], y[ind], color=c, linewidth=3, label = 'rr14-steep' if label == True else None)
        ind = np.where((yun < -1) & (zin == z))
        ax.plot(x[ind], yun[ind], color=c, linewidth=2, linestyle='dotted', label = 'unattenuated' if label == True else None)

    x, y, yun, ydi, ym, yd, zin = np.loadtxt("../data/Models/SharkVariations/uv_lf_L24HBTbestparams_eagle-rr14_L800_allz.txt", unpack = True, usecols = [0,1,2,3,4,5,6])

    if(all_channels == False):
        ind = np.where((y < -1) & (zin == z))
        ax.plot(x[ind], y[ind], color=c, linewidth=3, linestyle='dashed', label = 'rr14' if label == True else None)

plt = common.load_matplotlib()
outdir = '/scratch/pawsey0119/clagos/SHARK_Out/Plots/L800/Sharkv2-Lagos23-bestparams/eagle-rr14-steep/'
obsdir = '/software/projects/pawsey0119/clagos/shark/data/'

fig = plt.figure(figsize=(5,4))
xmin, xmax, ymin, ymax = -25, -13, -7, -1
xtitplot="$\\rm 1500\\AA\, mag\, (AB)$"
ytitplot="$\\rm log_{10}(\Phi/{\\rm dex^{-1}} {\\rm Mpc}^{-3})$"


idx = (2, 3, 4, 5, 6, 7)
zs  = (2, 3, 4, 5, 6, 7)
zs_label = ['6', '8', '10', '13', '15', '17']
zs_in = [6, 8, 10, 13, 15, 17]

cols = ('Navy','Aquamarine', 'Green','Gold','red','DarkRed') 

ax = fig.add_subplot(111)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))

for idx, z, c, zl in zip(idx, zs_in, cols, zs_label):
    plot_uv_lf_l800_variants(ax, zl, z, label = False, color=c)
ax.text(-24, -3, 'z=6', color=cols[0])
ax.text(-24, -3.4, 'z=8', color=cols[1])
ax.text(-24, -3.8, 'z=10', color=cols[2])
ax.text(-24, -4.2, 'z=13', color=cols[3])
ax.text(-24, -4.6, 'z=15', color=cols[4])
ax.text(-24, -5, 'z=17', color=cols[5])
ax.plot([-24,-22.8],[-1.5,-1.5],linestyle='solid', color='black')
ax.text(-22.7,-1.5,'fiducial')
ax.plot([-24,-22.8],[-2,-2],linestyle='dashed', color='black')
ax.text(-22.7,-2,'no disk ins.')
ax.plot([-24,-22.8],[-2.5,-2.5],linestyle='dotted', color='black')
ax.text(-22.7,-2.5,'$\\eta_{\\rm burst}=1$')



common.prepare_legend(ax, cols)

plt.tight_layout()
common.savefig(outdir, fig, "UV_luminosity_function_evolution_zGT6_nodiskins.pdf")

#dust scaling variants
fig = plt.figure(figsize=(5,4))
xmin, xmax, ymin, ymax = -25, -14, -7, -1
xtitplot="$\\rm 1500\\AA\, mag\, (AB)$"
ytitplot="$\\rm log_{10}(\Phi/{\\rm dex^{-1}} {\\rm Mpc}^{-3})$"


idx = (2, 3, 4, 5, 6, 7)
zs  = (2, 3, 4, 5, 6, 7)
zs_label = ['6', '8', '10', '13', '15', '17']
zs_in = [6, 8, 10, 13, 15, 17]

cols = ('Navy','Aquamarine', 'Green','Gold','red','DarkRed') 

ax = fig.add_subplot(111)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))

for idx, z, c, zl in zip(idx, zs_in, cols, zs_label):
    plot_uv_lf_l800_dust_variants(ax, zl, z, label = False, color=c)
ax.text(-24.5, -3, 'z=6', color=cols[0])
ax.text(-24.5, -3.4, 'z=8', color=cols[1])
ax.text(-24.5, -3.8, 'z=10', color=cols[2])
ax.text(-24.5, -4.2, 'z=13', color=cols[3])
ax.text(-24.5, -4.6, 'z=15', color=cols[4])
ax.text(-24.5, -5, 'z=17', color=cols[5])
ax.plot([-24.5,-23.5],[-2,-2],linestyle='solid', color='black')
ax.text(-23.3,-2,'rr14-steep (L800)')
ax.plot([-24.5,-23.5],[-2.5,-2.5],linestyle='dashed', color='black')
ax.text(-23.3,-2.5,'rr14 (L800)')
ax.plot([-24.5,-23.5],[-1.5,-1.5],linestyle='dotted', color='black')
ax.text(-23.3,-1.5,'unattenuated(L800)')


common.prepare_legend(ax, cols)

plt.tight_layout()
common.savefig(outdir, fig, "UV_luminosity_function_evolution_zGT6_dust_scaling.pdf")

fig = plt.figure(figsize=(12,8))
xmin, xmax, ymin, ymax = -25, -11, -7, -0.5


xtit="$\\rm 1500\\AA\, mag\, (AB)$"
ytit="$\\rm log_{10}(\Phi/{\\rm dex^{-1}} {\\rm Mpc}^{-3})$"

xleg = xmin + 0.2 * (xmax-xmin)
yleg = ymax - 0.1 * (ymax-ymin)

labels= ('z=3', 'z=4', 'z=6', 'z=8', 'z=10', 'z=13', 'z=15,17')

corrm_obs = -5.0*np.log10(0.67/0.7)
corry_obs = 3.0*np.log10(0.67/0.7)

subplots = (231, 232, 233, 234, 235)
idx = (2, 3, 4, 5, 6)
zs  = (2, 3, 4, 5, 6)

for subplot, idx, z in zip(subplots, idx, zs):

    ax = fig.add_subplot(subplot)
    if (idx > 3):
        xtitplot = xtit
    else:
        xtitplot = ' '
    if ((idx == 2) | (idx == 5)):
        ytitplot = ytit
    else:
        ytitplot = ' '

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
    ax.text(xleg,yleg, labels[idx])

    if(idx == 2):
       file = obsdir+'/lf/lf1600_z4-10_Bouwens2015.data'
       lmB15z4,p4,dp4,lmB15z6,p6,dp6,lmB15z8,p8,dp8 = np.loadtxt(file,usecols=[0, 1, 2, 6, 7, 8, 12, 13, 14],unpack=True)

       file = obsdir+'/lf/lf1500_z4-8_Finkelstein2015.data'
       lmF15,pF4,dpuF4,dpdF4,pF6,dpuF6,dpdF6,pF8,dpuF8,dpdF8 = np.loadtxt(file,usecols=[0,1, 2, 3, 7, 8, 9, 13, 14, 15],unpack=True)

       yobs = np.log10(p6)
       ydn  = np.log10(p6-dp6)
       yup  = np.log10(p6+dp6)
       ax.errorbar(lmB15z6+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^', label="Bouwens+2015")

       yobs = np.log10(pF6*1e-3)
       ydn  = np.log10(pF6*1e-3-dpdF6*1e-3)
       yup  = np.log10(pF6*1e-3+dpuF6*1e-3)
       ax.errorbar(lmF15+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='v', label="Finkelstein+2015")


       file = obsdir+'/lf/lf1500_z7_atek2026.dat'
       lma26,pa26,dpa26 = np.loadtxt(file,usecols=[0, 1, 2],unpack=True)
       ax.errorbar(lma26+corrm_obs, np.log10(pa26)+corry_obs, yerr=[np.log10(pa26)-np.log10(pa26-dpa26), np.log10(pa26+dpa26)-np.log10(pa26)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='P', label="Atek+2026")

       plot_uv_lf_micro(ax, 6)

       plot_uv_lf_hisurfs(ax, 6)
       plot_uv_lf_l800(ax, 6)
       plot_uv_lf_medi(ax, 6)

    if(idx == 3):
       yobs = np.log10(p8)
       ydn  = np.log10(p8-dp8)
       yup  = np.log10(p8+dp8)
       ax.errorbar(lmB15z8+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^')

       yobs = np.log10(pF8*1e-3)
       ydn  = np.log10(pF8*1e-3-dpdF8*1e-3)
       yup  = np.log10(pF8*1e-3+dpuF8*1e-3)
       ax.errorbar(lmF15+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='v')

       file = obsdir+'/lf/lf1500_z8_adams23.data'
       lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
       ax.errorbar(lm+corrm_obs, np.log10(p * 1e-5)+corry_obs, yerr=[np.log10(p)-np.log10(dp),np.log10(p + dp) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D', label='Adams+2024')

       ax.errorbar(lma26+corrm_obs, np.log10(pa26)+corry_obs, yerr=[np.log10(pa26)-np.log10(pa26-dpa26), np.log10(pa26+dpa26)-np.log10(pa26)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='P')

       plot_uv_lf_micro(ax, 8, label = False)
       plot_uv_lf_hisurfs(ax, 8, label = False)
       plot_uv_lf_l800(ax, 8, label = False)
       plot_uv_lf_medi(ax, 8, label = False)

    if(idx == 4):
       file = obsdir+'/lf/lf1500_z10_oesch2018.data'
       lm,p,dpu,dpd = np.loadtxt(file,usecols=[0,1, 2, 3],unpack=True)
       ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[p-dpu,dpd-p], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*', label='Oesch+2018')

       file = obsdir+'/lf/lf1500_z10_adams23.data'
       lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
       ax.errorbar(lm+corrm_obs, np.log10(p * 1e-5)+corry_obs, yerr=[np.log10(p)-np.log10(dp),np.log10(p + dp) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D')

       file = obsdir+'/lf/lf1500_z10_weibel25.data'
       lm,p,dpup,dpdn,flag = np.loadtxt(file,usecols=[0,1,2,3,4],unpack=True)
       ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[dpup, dpdn], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o', label='Weibel+2026')
       
       file = obsdir+'/lf/lf1500_z10_whitler25.data'
       lm,p,dpdn,dpup = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
       #ax.errorbar(lm+corrm_obs, np.log10(p)+corry_obs, yerr=[np.log10(p) - np.log10(p-dpdn), np.log10(p + dpup) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s', label='Whitler+25')
       plot_uv_lf_micro(ax, 10, label = False)
       plot_uv_lf_hisurfs(ax, 10, label = False)
       plot_uv_lf_l800(ax, 10, label = False)
       plot_uv_lf_medi(ax, 10, label = False)

    if(idx == 5):
       file = obsdir+'/lf/lf1500_z13_weibel25.data'
       lm,p,dpup,dpdn,flag = np.loadtxt(file,usecols=[0,1,2,3,4],unpack=True)
       ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[dpup, dpdn], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')
       file = obsdir+'/lf/lf1500_z12p8_whitler25.data'
       lm,p,dpdn,dpup = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
       #ax.errorbar(lm+corrm_obs, np.log10(p)+corry_obs, yerr=[np.log10(p) - np.log10(p-dpdn), np.log10(p + dpup) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')
       plot_uv_lf_micro(ax, 12.56, label = False)
       plot_uv_lf_hisurfs(ax, 12.62, label = False)
       plot_uv_lf_l800(ax, 13, label = False)
       plot_uv_lf_medi(ax, 13, label = False)

    if(idx == 6):
       file = obsdir+'/lf/lf1500_z17_weibel25.data'
       lm,p= np.loadtxt(file,usecols=[0,1],unpack=True)
       ax.plot(lm+corrm_obs, p+corry_obs, ls='None',marker='o', fillstyle ='none', color='grey')
       for x,y in zip(lm+corrm_obs, p+corry_obs):
           ax.arrow(x,y,0,-0.75, head_width=0.05, head_length=0.1, color='grey')
       file = obsdir+'/lf/lf1500_z17_perez-gonzalez25.data'
       lm,p,dn,du,flag= np.loadtxt(file,usecols=[0,1,2,3,4],unpack=True)
       flagin = np.where(flag == 0)
       ax.errorbar(lm[flagin]+corrm_obs+0.15, p[flagin]+corry_obs, yerr=[dn[flagin], du[flagin]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s', label = 'Perez-Gonzalez+2025')

       ax.plot(lm[0]+corrm_obs, p[0]+corry_obs, ls='None',marker='s', fillstyle ='none', color='grey')
       ax.arrow(lm[0]+corrm_obs, p[0]+corry_obs,0,-0.75, head_width=0.05, head_length=0.1, color='grey')

       file = obsdir+'/lf/lf1500_z15_whitler25.data'
       lm,p,dpdn,dpup = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
       #ax.errorbar(lm+corrm_obs, np.log10(p)+corry_obs, yerr=[np.log10(p) - np.log10(p-dpdn), np.log10(p + dpup) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')
       plot_uv_lf_micro(ax, 15, label = False)
       plot_uv_lf_micro(ax, 17, label = False)
       plot_uv_lf_hisurfs(ax, 15, label = False)
       plot_uv_lf_hisurfs(ax, 17, label = False)
       plot_uv_lf_l800(ax, 15, label = False)
       plot_uv_lf_l800(ax, 17, label = False)
       plot_uv_lf_medi(ax, 15, label = False)
       plot_uv_lf_medi(ax, 17, label = False)

    if(idx == 2):
        common.prepare_legend(ax, ['Purple','Purple','Orchid', 'Orchid', 'Olive','Olive','k','k','grey','grey','grey'],  bbox_to_anchor=(2.205, -1.22))
    if idx == 3:
        common.prepare_legend(ax, ['grey','grey','grey'], bbox_to_anchor=(1.01, -1.31))
    if idx == 4:
        common.prepare_legend(ax, ['grey','grey','grey'], bbox_to_anchor=(0.51, -1.13))
    if idx == 6:
        common.prepare_legend(ax, ['grey','grey','grey'], bbox_to_anchor=(1.71, -0.01))


    #if ((idx ==4) or (idx == 3)):
    #    common.prepare_legend(ax, ['grey','grey','grey'], loc=0)

#plt.tight_layout()
common.savefig(outdir, fig, "UV_luminosity_function_evolution_zGT6.pdf")

fig = plt.figure(figsize=(12,8))
xmin, xmax, ymin, ymax = -25, -13, 0, 1

subplots = (231, 232, 233, 234, 235)
idx = (2, 3, 4, 5, 6)
zl800  = (6, 8, 10, 12.7, 15)
zhisurfs  = (6, 8, 10, 12.62, 15)
zmicrosurfs  = (6, 8, 10, 12.56, 15)

ytit="Fractional channel contribution"
labels= ('z=3', 'z=4', 'z=6', 'z=8', 'z=10', 'z=13', 'z=15')

p=0
for subplot, idx in zip(subplots, idx):

    ax = fig.add_subplot(subplot)
    if (idx > 3):
        xtitplot = xtit
    else:
        xtitplot = ' '
    if ((idx == 2) | (idx == 5)):
        ytitplot = ytit
    else:
        ytitplot = ' '

    common.prepare_ax(ax, -23, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 0.1, 0.1))
    ax.text(-15,0.9, labels[idx])

    if(idx == 2):
       plot_uv_contributions_l800(ax, zl800[p], label = True)
       plot_uv_contributions_hisurfs(ax, zhisurfs[p], label = True)
       plot_uv_contributions_microsurfs(ax, zmicrosurfs[p], label = True)

    else:
       plot_uv_contributions_l800(ax, zl800[p], label = False)
       plot_uv_contributions_hisurfs(ax, zhisurfs[p], label = False)
       plot_uv_contributions_microsurfs(ax, zmicrosurfs[p], label = False)

    p = p + 1

    if(idx == 2):
        common.prepare_legend(ax, ['Navy','DarkSeaGreen','Firebrick', 'Navy','Navy','DarkSeaGreen','Firebrick', 'grey','grey','grey'],  bbox_to_anchor=(2.305, -0.8))
    #if(idx == 6):
    #    ax.text(-10,0.75, "Shark v2.0 (L800)", fontsize=12)

#plt.tight_layout()
common.savefig(outdir, fig, "UV_luminosity_function_channels_contribution_zGT6.pdf")


fig = plt.figure(figsize=(5,7))
ymin, ymax, xmin, xmax = 5, 10, -13.7, -23

xtit="$\\rm1500\\AA\, mag\, (AB)$"
ytit="$\\rm log_{10}(M_{\\star}/M_{\\odot})$"

labels= ('z=6', 'z=8', 'z=10', 'z=13', 'z=15', 'z=17')
zl800  = (6, 8, 10, 12.7, 15, 17)
zhisurfs  = (6, 8, 10, 12.62, 15, 17)
zmicrosurfs  = (6, 8, 10, 13, 15, 17)
uvlim = [-14, -14.5, -16, -16.5, -17, -17]

ax = fig.add_subplot(211)

common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

cols = ('Navy','Aquamarine', 'Green','Gold','red','DarkRed') 

for i in range(0,len(labels)):

    #plot_uv_mass_microsurfs(ax, zmicrosurfs[i], color = cols[i], stellar_mass = True, scatter=False, label = True, labelin=labels[i])
    plot_uv_mass_l800(ax, zl800[i], color = cols[i], stellar_mass = True, scatter=False, label = True, labelin=labels[i], uvlim = uvlim[i])
    plot_uv_mass_l800_nodiskins(ax, zl800[i], color = cols[i], stellar_mass = True, scatter=False, label = True, labelin=labels[i], uvlim = uvlim[i])
    #plot_uv_mass_hisurfs(ax, zhisurfs[i], color = cols[i], stellar_mass = True, scatter=False, label = True, labelin=labels[i])
    ax.text(-22, 9.1 - 0.4 * i, labels[i], color=cols[i])

ax.text(-19, 6.2, 'Shark v2.0')
ax.plot([-19,-19.7],[5.8, 5.8], linestyle='solid', linewidth=2, color='k')
ax.text(-19.8,5.8, 'fiducial (L800)')
ax.plot([-19,-19.7],[5.4,5.4], linestyle='dashed', linewidth=1, color='k')
ax.text(-19.8,5.4, 'no disk ins. (L800)')
#ax.plot([8.2,8.7],[-12.4,-12.4], linestyle='solid', linewidth=1, color='k')
#ax.text(8.8,-12.4, 'fiducial (micro)')
common.prepare_legend(ax, cols, loc=0)

ax = fig.add_subplot(212)

ytit="$\\rm \\sigma(log_{10}(M_{\\star}/M_{\\odot}))$"

ymin, ymax = -0.6,0.6

common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 0.1, 0.1))

for i in range(0,len(labels)):

    #plot_uv_mass_microsurfs(ax, zmicrosurfs[i], color = cols[i], stellar_mass = True, scatter=True, label = False)
    plot_uv_mass_l800(ax, zl800[i], color = cols[i], stellar_mass = True, scatter=True, label = False, uvlim = uvlim[i])
    plot_uv_mass_l800_nodiskins(ax, zl800[i], color = cols[i], stellar_mass = True, scatter=True, label = False, uvlim = uvlim[i])
    #plot_uv_mass_hisurfs(ax, zhisurfs[i], color = cols[i], stellar_mass = True, scatter=True, label = False)

ax.plot([xmin, xmax], [0,0], linestyle='dotted', color='k')

plt.tight_layout()
common.savefig(outdir, fig, "UV_StellarMass_zGT6.pdf")

fig = plt.figure(figsize=(5,7))
ymin, ymax, xmin, xmax = 9.4, 12, -13, -23

xtit="$\\rm 1500\\AA\, mag\, (AB)$"
ytit="$\\rm log_{10}(M_{\\rm halo}/M_{\\odot})$"
#idxs = (2, 3, 4, 5)
#zs  = (2, 3, 4, 5)

ax = fig.add_subplot(211)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 0.5, 0.5))

for i in range(0,len(labels)):

    #Predicted LF
    #plot_uv_mass_microsurfs(ax, zmicrosurfs[i], color = cols[i], stellar_mass = False, scatter=False, label = True, labelin = labels[i])
    plot_uv_mass_l800(ax, zl800[i], color = cols[i], stellar_mass = False, scatter=False, label = False, uvlim = uvlim[i])
    plot_uv_mass_l800_nodiskins(ax, zl800[i], color = cols[i], stellar_mass = False, scatter=False, label = False, uvlim = uvlim[i])

    #plot_uv_mass_hisurfs(ax, zhisurfs[i], color = cols[i], stellar_mass = False, scatter=False, label = False)

common.prepare_legend(ax, cols, loc=0)

ax = fig.add_subplot(212)
ytit="$\\rm \\sigma(log_{10}(M_{\\rm halo}/M_{\\odot}))$"

ymin, ymax = -0.6, 0.6
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 0.1, 0.1))

for i in range(0,len(labels)):

    #Predicted LF
    #plot_uv_mass_microsurfs(ax, zmicrosurfs[i], color = cols[i], stellar_mass = False, scatter=True, label = False)
    plot_uv_mass_l800(ax, zl800[i], color = cols[i], stellar_mass = False, scatter=True, label = False, uvlim = uvlim[i])
    plot_uv_mass_l800_nodiskins(ax, zl800[i], color = cols[i], stellar_mass = False, scatter=True, label = False, uvlim = uvlim[i])

    #plot_uv_mass_hisurfs(ax, zhisurfs[i], color = cols[i], stellar_mass = False, scatter=True, label = False)

ax.plot([xmin, xmax], [0,0], linestyle='dotted', color='k')

plt.tight_layout()
common.savefig(outdir, fig, "UV_HaloMass_zGT6.pdf")



