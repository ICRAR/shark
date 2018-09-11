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
"""Size plots"""

import functools

import numpy as np

import common
import utilities_statistics as us

# Initialize arguments
zlist = (0, 0.5, 1, 2)

##################################
#Constants
RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun

mlow = 6.5
mupp = 12.5
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

dmobs = 0.4
mbins_obs = np.arange(mlow,mupp,dmobs)
xmf_obs = mbins_obs + dmobs/2.0

def prepare_data(hdf5_data, index, rcomb, disk_size, bulge_size, bulge_size_mergers, bulge_size_diskins, BH,
                 disk_size_sat, disk_size_cen, BT_fractions, BT_fractions_nodiskins, bulge_vel, 
                 disk_vel, BT_fractions_centrals, BT_fractions_satellites):

    (h0, _, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     mBH, rdisk, rbulge, typeg, specific_angular_momentum_disk_star, specific_angular_momentum_bulge_star, 
     specific_angular_momentum_disk_gas, specific_angular_momentum_bulge_gas, specific_angular_momentum_disk_gas_atom, 
     specific_angular_momentum_disk_gas_mol, lambda_sub, mvir_s) = hdf5_data

    mbulge_mergers = mburst_mergers + mstars_bulge_mergers_assembly
    zero_bulge = np.where(rbulge <= 0)
    if(len(rbulge) == len(rbulge[zero_bulge])):
            #case where there is zero bulge build up.
            rbulge[zero_bulge] = 1e-10
            specific_angular_momentum_bulge_star[zero_bulge] = 1.0
            mbulge[zero_bulge] = 10.0

    bin_it = functools.partial(us.wmedians, xbins=xmf)

    vdisk = specific_angular_momentum_disk_star / rdisk / 2.0 #in km/s
    vbulge = specific_angular_momentum_bulge_star / rbulge / 2.0 #in km/s
   
    ind = np.where(mdisk+mbulge > 0)
    rcomb[index,:] = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                            y=np.log10((mdisk[ind]*rdisk[ind]  + mbulge[ind]*rbulge[ind])*MpcToKpc / (mdisk[ind]+mbulge[ind])))
    BT_fractions[index] = us.fractional_contribution(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),y=mbulge[ind]/(mdisk[ind]+mbulge[ind]), xbins=xmf)

    BT_fractions_nodiskins[index] = us.fractional_contribution(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
		                    y=(mbulge_mergers[ind])/(mdisk[ind]+mbulge[ind]), xbins=xmf)

    ind = np.where((mdisk+mbulge > 0) & (typeg == 0))
    BT_fractions_centrals[index] = us.fractional_contribution(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),y=mbulge[ind]/(mdisk[ind]+mbulge[ind]), xbins=xmf)
    ind = np.where((mdisk+mbulge > 0) & (typeg > 0))
    BT_fractions_satellites[index] = us.fractional_contribution(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),y=mbulge[ind]/(mdisk[ind]+mbulge[ind]), xbins=xmf)

    ind = np.where((mdisk > 0)  & (mdisk/(mdisk+mbulge) > 0.5))
    disk_size[index,:] = bin_it(x=np.log10(mdisk[ind]) - np.log10(float(h0)),
                                y=np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)))

    ind = np.where((mdisk > 0) & (typeg == 0) & (mbulge/mdisk < 0.5))
    disk_vel[index,:] = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                y=np.log10(vdisk[ind]))

    ind = np.where((mdisk > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) > 0.5))
    disk_size_cen[index,:]  = bin_it(x=np.log10(mdisk[ind]) - np.log10(float(h0)),
                                    y=np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)))

    ind = np.where((mdisk > 0) & (typeg > 0) & (mdisk/(mdisk+mbulge) > 0.5))
    disk_size_sat[index,:] = bin_it(x=np.log10(mdisk[ind]) - np.log10(float(h0)),
                                    y=np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)))

    ind = np.where((mbulge > 0) & (mbulge/(mbulge+mdisk) > 0.5) & (rbulge > 1e-6))
    bulge_size[index,:] = bin_it(x=np.log10(mbulge[ind]) - np.log10(float(h0)),
                                 y=np.log10(rbulge[ind]*MpcToKpc) - np.log10(float(h0)))

    ind = np.where((mbulge > 0) & (mbulge/(mbulge+mdisk) > 0.5) & (rbulge > 1e-6) & (mbulge_mergers/mbulge > 0.5))
    bulge_size_mergers[index,:] = bin_it(x=np.log10(mbulge[ind]) - np.log10(float(h0)),
                                 y=np.log10(rbulge[ind]*MpcToKpc) - np.log10(float(h0)))

    ind = np.where((mbulge > 0) & (mbulge/(mbulge+mdisk) > 0.5) & (rbulge > 1e-6) & ((mbulge - mbulge_mergers) / mbulge > 0.5))
    bulge_size_diskins[index,:] = bin_it(x=np.log10(mbulge[ind]) - np.log10(float(h0)),
                                 y=np.log10(rbulge[ind]*MpcToKpc) - np.log10(float(h0)))

    ind = np.where(mbulge > 0)
    BH[index,:] = bin_it(x=np.log10(mbulge[ind]) - np.log10(float(h0)),
                    y=np.log10(mBH[ind]) - np.log10(float(h0)))
    
    ind = np.where((mbulge > 0) & (mbulge/mdisk > 0.5))
    bulge_vel[index,:] = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                    y=np.log10(vbulge[ind]))
    

def plot_sizes(plt, outdir, obsdir, disk_size_cen, disk_size_sat, bulge_size, bulge_size_mergers, bulge_size_diskins):

    fig = plt.figure(figsize=(5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star,disk}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm r_{\\star,disk}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.1, 2
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # LTG ##################################
    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(disk_size_cen[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = disk_size_cen[0,0,ind]
    errdn = disk_size_cen[0,1,ind]
    errup = disk_size_cen[0,2,ind]
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="Shark centrals")

    #Predicted size-mass for disks in disk=dominated galaxies satellites
    ind = np.where(disk_size_sat[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = disk_size_sat[0,0,ind]
    errdn = disk_size_sat[0,1,ind]
    errup = disk_size_sat[0,2,ind]
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'r', mec='r',marker='v',markersize=5, label="Shark satellites")

    #Lange et al. (2016)
    a = 5.56
    aerr = 1.745
    b = 0.274
    ind = np.where(xmf < 10.3)
    rL16 = np.log10(a*pow((pow(10.0,xmf)/1e10),b))

    ax.plot(xmf[ind],rL16[ind],'b', linestyle='solid', label ='L16 disks')


    m,r = common.load_observation(obsdir, 'SizesAndAM/rdisk_L16.dat', [0,1])
    ax.plot(m[0:36], r[0:36], linestyle='dotted',color='b',label="50th, 68th, 90th")
    ax.plot(m[38:83], r[38:83], linestyle='dotted',color='b')
    ax.plot(m[85:128], r[85:129], linestyle='dotted',color='b')

    common.prepare_legend(ax, ['b','b','k','r'], loc=2)

    # ETGs ##################################
    xtit = "$\\rm log_{10} (\\rm M_{\\star,bulge}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm r_{\\star,bulge}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.35, 2
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(212)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #Predicted size-mass for bulges in bulge-dominated systems
    ind = np.where(bulge_size[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = bulge_size[0,0,ind]
        errdn = bulge_size[0,1,ind]
        errup = bulge_size[0,2,ind]
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="Shark bulges")
        #for i in zip(bulge_size[0,0,:]):
        #    print i

    ind = np.where((bulge_size_diskins[0,0,:] != 0) & (xmf > 10.2))
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = bulge_size_diskins[0,0,ind]
        err   = bulge_size[0,1,ind]
        err[:] = 0
        ax.errorbar(xplot,yplot[0],yerr=[err[0],err[0]], ls='None', mfc='None', ecolor = 'Orange', mec='Orange',marker='s', markersize=4, label="disk instability driven")

    ind = np.where((bulge_size_mergers[0,0,:] != 0) & (xmf > 10.2))
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = bulge_size_mergers[0,0,ind]
        err   = bulge_size[0,1,ind]
        err[:] = 0
        ax.errorbar(xplot,yplot[0],yerr=[err[0],err[0]], ls='None', mfc='None', ecolor = 'DarkCyan', mec='DarkCyan',marker='D',  markersize=4, label="merger driven")

    rb_nodissipation = common.load_observation(obsdir, 'Models/SharkVariations/SizeBulges_OtherModels.dat', [0])
    ind = np.where(rb_nodissipation != 0)
    xplot = xmf[ind]
    yplot = rb_nodissipation[ind]
    err   = xmf[ind]
    err[:] = 0
    ax.errorbar(xplot,yplot,yerr=[err,err], ls='None', mfc='None', ecolor = 'LightSlateGray', mec='LightSlateGray',marker='x',  markersize=6, label="no dissipation")

    #Lange et al. (2016)
    a = 2.319
    aerr = 1.186
    b = 0.19565217391
    a2 = 1.390
    aerr2 = 0.9
    b2 = 0.624

    ind = np.where(xmf <= 10.3)

    rL16 = np.log10(a*pow((pow(10.0,xmf[ind])/1e10),b))
    ax.plot(xmf[ind],rL16,'r', linestyle='solid', label ='L16 $M_{\\star}<2\\times 10^{10}\\rm M_{\odot}$')

    ind = np.where((xmf >= 10.3) & (xmf < 11.5))
    rL16_2 = np.log10(a2*pow((pow(10.0,xmf[ind])/1e10),b2))

    ax.plot(xmf[ind],rL16_2,'m', linestyle='solid', label ='L16 $M_{\\star}>2\\times 10^{10}\\rm M_{\odot}$')
    
    m,r = common.load_observation(obsdir, 'SizesAndAM/rbulge_L16.dat', [0,1])
    ax.plot(m[0:39], r[0:39], linestyle='dotted',color='r')
    ax.plot(m[41:76], r[41:76], linestyle='dotted',color='r')
    ax.plot(m[78:115], r[78:115], linestyle='dotted',color='r')

    common.prepare_legend(ax, ['r','m','k','Orange','DarkCyan','LightSlateGray'], loc=2)
    common.savefig(outdir, fig, 'sizes.pdf')


def plot_velocities(plt, outdir, disk_vel, bulge_vel):

    fig = plt.figure(figsize=(5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm v_{\\rm 50, disk}/km s^{-1})$"
    xmin, xmax, ymin, ymax = 8, 12, 1, 3
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # LTG ##################################
    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'z=0')

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(disk_vel[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = disk_vel[0,0,ind]
    errdn = disk_vel[0,1,ind]
    errup = disk_vel[0,2,ind]
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="Shark centrals B/T<0.5")

    common.prepare_legend(ax, ['k'], loc=2)

    # ETGs ##################################
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm v_{\\rm 50, bulge}/km s^{-1})$"
    xmin, xmax, ymin, ymax = 8, 12, 1, 3
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(212)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'z=0')

    #Predicted size-mass for bulges in bulge-dominated systems
    ind = np.where(bulge_vel[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = bulge_vel[0,0,ind]
        errdn = bulge_vel[0,1,ind]
        errup = bulge_vel[0,2,ind]
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="Shark bulges B/T > 0.5")

        common.prepare_legend(ax, ['k'], loc=2)

    common.savefig(outdir, fig, 'velocities.pdf')
 
def plot_sizes_combined(plt, outdir, rcomb):

    fig = plt.figure(figsize=(5,4.5))

    # Total ##################################
    xtit="$\\rm log_{10} (\\rm M_{\\rm stars, total}/M_{\odot})$"
    ytit="$\\rm log_{10} (\\rm r_{\\rm 50, comb}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.5, 2

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    #Predicted size-mass for disks
    ind = np.where(rcomb[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = rcomb[0,0,ind]
    errdn = rcomb[0,1,ind]
    errup = rcomb[0,2,ind]
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="Shark disk+bulge combined")

    common.prepare_legend(ax, ['k','k','k'], loc=2)
    common.savefig(outdir, fig, 'sizes_combined.pdf')


def plot_bulge_BH(plt, outdir, obsdir, BH):

    fig = plt.figure(figsize=(5,4.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm bulge}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm M_{\\rm BH}/M_{\odot})$"

    xmin, xmax, ymin, ymax = 8, 12, 5, 11
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'z=0')

    #Predicted SMHM
    ind = np.where(BH[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = BH[0,0,ind]
        errdn = BH[0,1,ind]
        errup = BH[0,2,ind]
        ax.plot(xplot,yplot[0],color='k',label="Shark")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='grey', interpolate=True)


    MBH_othermodels = common.load_observation(obsdir, 'Models/SharkVariations/BHBulgeRelation_OtherModels.dat', [0])
    MBH_f_smbh0p008   = MBH_othermodels[0:29]
    MBH_f_smbh0p00008 = MBH_othermodels[30:60]
    ind = np.where(MBH_f_smbh0p008 != 0)
    xplot = xmf[ind]
    yplot = MBH_f_smbh0p008[ind]
    ax.plot(xplot,yplot,color='Goldenrod',linestyle='dashed',label='$f_{\\rm smbh}=8 \\times 10^{-3}$')
    ind = np.where(MBH_f_smbh0p00008 != 0)
    xplot = xmf[ind]
    yplot = MBH_f_smbh0p00008[ind]
    ax.plot(xplot,yplot,color='Orange',linestyle='dotted',label='$f_{\\rm smbh}=8 \\times 10^{-5}$')

    #BH-bulge relation
    mBH_M13, errup_M13, errdn_M13, mBH_power, mbulge_M13 = common.load_observation(obsdir, 'BHs/MBH_sigma_Mbulge_McConnelMa2013.dat', [0,1,2,3,7])

    ind = np.where((mBH_M13 > 0) & (mbulge_M13 > 0))
    xobs = np.log10(mbulge_M13[ind])
    yobs = np.log10(mBH_M13[ind] * pow(10.0,mBH_power[ind]))
    lerr = np.log10((mBH_M13[ind] - errdn_M13[ind]) * pow(10.0,mBH_power[ind]))
    herr = np.log10((mBH_M13[ind] + errup_M13[ind]) * pow(10.0,mBH_power[ind]))
    ax.errorbar(xobs, yobs, yerr=[yobs-lerr,herr-yobs], ls='None', mfc='None', ecolor = 'r', mec='r',marker='^',label="McConnell & Ma 2013")

    #BH-bulge relation
    mBH_H04, errup_H04, errdn_H04, mbulge_H04 = common.load_observation(obsdir, 'BHs/MBH_sigma_Mbulge_HaeringRix2004.dat', [0,1,2,4])

    xobs = np.log10(mbulge_H04)

    yobs = xobs*0. - 999.
    indx = np.where( mBH_H04 > 0)
    yobs[indx] = np.log10(mBH_H04[indx])

    lerr = yobs*0. - 999.
    indx = np.where( (mBH_H04-errdn_H04) > 0)
    lerr[indx]  = np.log10(mBH_H04[indx] - errdn_H04[indx])

    herr = yobs*0. + 999.
    indx = np.where( (mBH_H04+errup_H04) > 0)
    herr[indx]  = np.log10(mBH_H04[indx] + errup_H04[indx])
    ax.errorbar(xobs, yobs, yerr=[yobs-lerr,herr-yobs], ls='None', mfc='None', ecolor = 'maroon', mec='maroon',marker='s',label="Haering+04")

    common.prepare_legend(ax, ['k','Goldenrod','Orange','r','maroon'], loc=2)
    common.savefig(outdir, fig, 'bulge-BH.pdf')


def plot_bt_fractions(plt, outdir, obsdir, BT_fractions, BT_fractions_nodiskins, BT_fractions_centrals, BT_fractions_satellites):

    fig = plt.figure(figsize=(5,4.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm f_{\\rm bulge}$"
    xmin, xmax, ymin, ymax = 8, 12, -0.1, 1.05

    # LTG ##################################
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #Predicted size-mass for disks
    ind = np.where(BT_fractions[0,:] >= 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = BT_fractions[0,ind]
        ax.plot(xplot,yplot[0],'k', label ='bulges built by all processes')

    ind = np.where(BT_fractions_nodiskins[0,:] >= 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = BT_fractions_nodiskins[0,ind]
        ax.plot(xplot,yplot[0],'r', linestyle = 'dashed', label ='only by mergers')

    #ind = np.where(BT_fractions_centrals[0,:] >= 0)
    #if(len(xmf[ind]) > 0):
    #    xplot = xmf[ind]
    #    yplot = BT_fractions_centrals[0,ind]
    #    ax.plot(xplot,yplot[0],'k', linestyle = 'dotted', label ='centrals')
    #ind = np.where(BT_fractions_satellites[0,:] >= 0)
    #if(len(xmf[ind]) > 0):
    #    xplot = xmf[ind]
    #    yplot = BT_fractions_satellites[0,ind]
    #    ax.plot(xplot,yplot[0],'k', linestyle = 'dashed', label ='satellites')

    BT_othermodels = common.load_observation(obsdir, 'Models/SharkVariations/BTFractions_OtherModels.dat', [0])
    BT_stable0   = BT_othermodels[0:29]
    BT_stable0p5 = BT_othermodels[30:60]
    BT_stable1   = BT_othermodels[91:120]

    ind = np.where(BT_stable0 >= 0)
    xplot = xmf[ind]
    yplot = BT_stable0[ind]
    ax.plot(xplot,yplot,color='Goldenrod',linestyle='dashdot',label='$\\epsilon_{\\rm disk}=0$')
    ind = np.where(BT_stable0p5 >= 0)
    xplot = xmf[ind]
    yplot = BT_stable0p5[ind]
    ax.plot(xplot,yplot,color='Orange',linestyle='dotted',label='$\\epsilon_{\\rm disk}=0.5$')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    mM16, fM16, errdnfM16, errupfM16 = common.load_observation(obsdir, 'Morph/Moffet16.dat', [0,1,2,3])
    errdnfM16 = np.abs(errdnfM16-fM16)
    errupfM16 = np.abs(errupfM16-fM16)
    ax.errorbar(mM16,fM16,yerr=[errdnfM16,errupfM16], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^',label="Moffett+16")

    common.prepare_legend(ax, ['k', 'r','Goldenrod', 'Orange','grey'], loc=2)
    common.savefig(outdir, fig, 'BTfractions.pdf')


def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'm_bh', 'rstar_disk', 'rstar_bulge', 'type', 
                           'specific_angular_momentum_disk_star', 'specific_angular_momentum_bulge_star',
                           'specific_angular_momentum_disk_gas', 'specific_angular_momentum_bulge_gas',
                           'specific_angular_momentum_disk_gas_atom', 'specific_angular_momentum_disk_gas_mol',
                           'lambda_subhalo', 'mvir_subhalo')}

    # Loop over redshift and subvolumes
    rcomb = np.zeros(shape = (len(zlist), 3, len(xmf)))
    disk_size = np.zeros(shape = (len(zlist), 3, len(xmf)))
    bulge_size = np.zeros(shape = (len(zlist), 3, len(xmf)))
    bulge_size_mergers = np.zeros(shape = (len(zlist), 3, len(xmf)))
    bulge_size_diskins = np.zeros(shape = (len(zlist), 3, len(xmf)))

    BH = np.zeros(shape = (len(zlist), 3, len(xmf)))
    disk_size_sat = np.zeros(shape = (len(zlist), 3, len(xmf)))
    disk_size_cen = np.zeros(shape = (len(zlist), 3, len(xmf)))
    BT_fractions = np.zeros(shape = (len(zlist), len(xmf)))
    BT_fractions_nodiskins = np.zeros(shape = (len(zlist), len(xmf)))
    BT_fractions_centrals = np.zeros(shape = (len(zlist), len(xmf)))
    BT_fractions_satellites = np.zeros(shape = (len(zlist), len(xmf)))
    disk_vel =  np.zeros(shape = (len(zlist), 3, len(xmf))) 
    bulge_vel =  np.zeros(shape = (len(zlist), 3, len(xmf)))
    
    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        prepare_data(hdf5_data, index, rcomb, disk_size, bulge_size, bulge_size_mergers, bulge_size_diskins, BH,
                     disk_size_sat, disk_size_cen, BT_fractions, BT_fractions_nodiskins, bulge_vel, disk_vel, 
                     BT_fractions_centrals, BT_fractions_satellites)

    plot_sizes(plt, outdir, obsdir, disk_size_cen, disk_size_sat, bulge_size, bulge_size_mergers, bulge_size_diskins)
    plot_velocities(plt, outdir, disk_vel, bulge_vel)
    plot_sizes_combined(plt, outdir, rcomb)
    plot_bulge_BH(plt, outdir, obsdir, BH)
    plot_bt_fractions(plt, outdir, obsdir, BT_fractions, BT_fractions_nodiskins, BT_fractions_centrals, BT_fractions_satellites)


if __name__ == '__main__':
    main(*common.parse_args())
