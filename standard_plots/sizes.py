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

import os
import common
import utilities_statistics as us

# Initialize arguments
zlist = (0, 0.5, 1, 2)

##################################
#Constants
RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun

mlow = 5
mupp = 12.5
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

dmobs = 0.4
mbins_obs = np.arange(mlow,mupp,dmobs)
xmf_obs = mbins_obs + dmobs/2.0

vlow = 1.0
vupp = 3.0
dv   = 0.1
vbins = np.arange(vlow,vupp,dv)
xv    = vbins + dv/2.0


def prepare_data(hdf5_data, index, rcomb, disk_size, bulge_size, bulge_size_mergers, bulge_size_diskins, BH,
                 disk_size_sat, disk_size_cen, BT_fractions, BT_fractions_nodiskins, bulge_vel, 
                 disk_vel, BT_fractions_centrals, BT_fractions_satellites, baryonic_TF, BHSM, xmf, xv, 
                 bs_error, BHSFR, BH_morpho, BHSM_morpho, age_stellar_mass):

    (h0, _, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     mBH, rdisk, rbulge, typeg, specific_angular_momentum_disk_star, specific_angular_momentum_bulge_star, 
     specific_angular_momentum_disk_gas, specific_angular_momentum_bulge_gas, specific_angular_momentum_disk_gas_atom, 
     specific_angular_momentum_disk_gas_mol, lambda_sub, mvir_s, mgas_disk, mgas_bulge, matom_disk, mmol_disk, matom_bulge, 
     mmol_bulge, mbh_acc_hh, mbh_acc_sb, age, sfr_disk, sfr_burst) = hdf5_data

    mstars_tot = (mdisk+mbulge)/h0
    #if index in (2, 3):
    #   for x, y, z, m in zip(mBH, mbh_acc_hh, mbh_acc_sb, mstars_tot):
    #       if x > 1e5 and m > 1e8:
    #          print(x/h0, y/h0/1e9, z/h0/1e9)

    mbulge_mergers = mburst_mergers + mstars_bulge_mergers_assembly
    zero_bulge = np.where(rbulge <= 0)
    if(len(rbulge) == len(rbulge[zero_bulge])):
            #case where there is zero bulge build up.
            rbulge[zero_bulge] = 1e-10
            specific_angular_momentum_bulge_star[zero_bulge] = 1.0
            mbulge[zero_bulge] = 10.0

    bin_it   = functools.partial(us.wmedians, xbins=xmf, nmin=8)
    bin_it_v = functools.partial(us.wmedians, xbins=xv)

    vdisk = specific_angular_momentum_disk_star / rdisk / 2.0  #in km/s
    vbulge = specific_angular_momentum_bulge_star / rbulge / 2.0 #in km/s
 
    ind = np.where(mdisk+mbulge > 0)
    age_stellar_mass[index,:] = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
            y=age[ind])

    ind = np.where(mdisk+mbulge > 0)
    rcomb[index,:] = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                            y=np.log10((mdisk[ind]*rdisk[ind]  + mbulge[ind]*rbulge[ind])*MpcToKpc / (mdisk[ind]+mbulge[ind]))- np.log10(float(h0)))

    bs_error[index] = us.bootstrap_error(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                            y=np.log10((mdisk[ind]*rdisk[ind]  + mbulge[ind]*rbulge[ind])*MpcToKpc / (mdisk[ind]+mbulge[ind])) - np.log10(float(h0)), xbins=xmf)

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

    ind = np.where((mdisk > 0) & ((mgas_disk+mgas_bulge)/(mdisk+mbulge) > 1))
    baryonic_TF[index,:] = bin_it_v(x=np.log10(vdisk[ind]), 
                                y=np.log10(mdisk[ind]+mbulge[ind]+matom_disk[ind]+mmol_disk[ind]+
                                matom_bulge[ind]+mmol_bulge[ind]) - np.log10(float(h0)))

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

    ind = np.where(mBH > 0)
    BH[index,:] = bin_it(x=np.log10(mbulge[ind]) - np.log10(float(h0)),
                    y=np.log10(mBH[ind]) - np.log10(float(h0)))
    ind = np.where((mBH > 0) & (mbulge/(mdisk+mbulge) > 0.5))
    BH_morpho[index,0,:] = bin_it(x=np.log10(mbulge[ind]) - np.log10(float(h0)),
                    y=np.log10(mBH[ind]) - np.log10(float(h0)))
    ind = np.where((mBH > 0) & (mbulge/(mdisk+mbulge) <= 0.5))
    BH_morpho[index,1,:] = bin_it(x=np.log10(mbulge[ind]) - np.log10(float(h0)),
                    y=np.log10(mBH[ind]) - np.log10(float(h0)))

    ind = np.where(mBH > 0)
    BHSM[index,:] = bin_it(x=np.log10(mbulge[ind] + mdisk[ind]) - np.log10(float(h0)),
                    y=np.log10(mBH[ind]) - np.log10(float(h0)))
   
    ind = np.where((mBH > 0) & (mbulge/(mdisk+mbulge) > 0.5))
    BHSM_morpho[index,0,:] = bin_it(x=np.log10(mbulge[ind] + mdisk[ind]) - np.log10(float(h0)),
                    y=np.log10(mBH[ind]) - np.log10(float(h0)))
    ind = np.where((mBH > 0) & (mbulge/(mdisk+mbulge) <= 0.5))
    BHSM_morpho[index,1,:] = bin_it(x=np.log10(mbulge[ind] + mdisk[ind]) - np.log10(float(h0)),
                    y=np.log10(mBH[ind]) - np.log10(float(h0)))

    ssfr = (sfr_disk + sfr_burst) / 1e9 / (mbulge + mdisk)
    #apply lower limit to ssfr
    ind = np.where(ssfr < 1e-14)
    ssfr[ind] = 2e-14

    ind = np.where((mBH > 0) & ((mbulge + mdisk)/h0 > 1e10) & (typeg <= 0))
    BHSFR[index,:] = bin_it(x=np.log10(mBH[ind]) - np.log10(float(h0)), 
                            y=np.log10(ssfr[ind]))

    ind = np.where((mbulge > 0) & (mbulge/mdisk > 0.5))
    bulge_vel[index,:] = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                    y=np.log10(vbulge[ind]))
    
    ind = np.where(mdisk >= 5e7)

    return(np.log10(mdisk[ind]) - np.log10(float(h0)), np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)), 13.7969 - age[ind])
 
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
    m = np.log10(m)
    r = np.log10(r)

    ax.plot(m[0:39], r[0:39], linestyle='dotted',color='b',label="50th, 68th, 90th")
    ax.plot(m[39:78], r[39:78], linestyle='dotted',color='b')
    ax.plot(m[78:len(r)], r[78:len(r)], linestyle='dotted',color='b')

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
    m = np.log10(m)
    r = np.log10(r)
    ax.plot(m[0:50], r[0:50], linestyle='dotted',color='r')
    ax.plot(m[50:90], r[50:90], linestyle='dotted',color='r')
    ax.plot(m[90:len(r)], r[90:len(r)], linestyle='dotted',color='r')

    common.prepare_legend(ax, ['r','m','k','Orange','DarkCyan','LightSlateGray'], loc=2)
    common.savefig(outdir, fig, 'sizes.pdf')


def plot_velocities(plt, outdir, disk_vel, bulge_vel, baryonic_TF):

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
 
    fig = plt.figure(figsize=(5,5))
    ytit = "$\\rm log_{10} (\\rm M_{\\rm bar}/M_{\odot})$"
    xtit = "$\\rm log_{10} (\\rm V_{\\rm max}/km\, s^{-1})$"
    xmin, xmax, ymin, ymax = 1.3, 2.5, 6, 13
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # LTG ##################################
    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    ind = np.where(xv < 2.477)
    ax.plot(xv,3.94*xv+1.79,linestyle='solid',color='b',label='McGaugh 2012')
    ax.fill_between(xv,3.94*xv+1.8,3.94*xv+1.79-0.26-0.25, facecolor='PaleTurquoise', alpha=0.5, interpolate=True)
    ax.fill_between(xv,3.94*xv+1.8,3.94*xv+1.79+0.22+0.25, facecolor='PaleTurquoise', alpha=0.5, interpolate=True)
    ax.fill_between(xv,3.94*xv+1.8,3.94*xv+1.79-0.26, facecolor='LightSeaGreen', alpha=0.5, interpolate=True)
    ax.fill_between(xv,3.94*xv+1.8,3.94*xv+1.79+0.26, facecolor='LightSeaGreen', alpha=0.5, interpolate=True)

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(baryonic_TF[0,0,:] != 0)
    xplot = xv[ind]
    yplot = baryonic_TF[0,0,ind]
    errdn = baryonic_TF[0,1,ind]
    errup = baryonic_TF[0,2,ind]
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="Shark gas-dominated")

    common.prepare_legend(ax, ['b','k'], loc=2)

    common.savefig(outdir, fig, 'baryon-TF.pdf')

def plot_sizes_combined(plt, outdir, obsdir, rcomb):

    lm, lr, count, bs_err = common.load_observation(obsdir, 'SizeMass/GAMA_H-band_dlogM_0.25_reff.txt', [0,1,2,3])

    fig = plt.figure(figsize=(5,4.5))

    # Total ##################################
    xtit="$\\rm log_{10} (\\rm M_{\\rm stars, total}/M_{\odot})$"
    ytit="$\\rm log_{10} (\\rm med(r_{50})/kpc)$"
    xmin, xmax, ymin, ymax = 8, 11.5, -0.1, 1

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    #Predicted size-mass for disks
    ind = np.where(rcomb[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = rcomb[0,0,ind]
    errdn = rcomb[0,1,ind]
    errup = rcomb[0,2,ind]
    #np.save(os.path.join(outdir,'sizemass.npy'), np.array([xplot, yplot[0]]))
    ax.plot(xplot,yplot[0], color = '#2A9D8F', 
            linewidth=2,label="Shark galaxies")
    
    # Add GAMA H-band observations with bootstrapped error
    ax.errorbar(lm, lr,yerr = [bs_err, bs_err], marker ='v',
             ls = 'none', mfc = 'None', markersize=5, 
             color = 'gray', label = 'Lange+2015')

    common.prepare_legend(ax, ['k','k','k'], loc=2)
    common.savefig(outdir, fig, 'sizes_combined.pdf')


def plot_bulge_BH(plt, outdir, obsdir, BH, BHSM, BHSFR, BH_morpho, BHSM_morpho):

    fig = plt.figure(figsize=(6,5.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm bulge}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm M_{\\rm BH}/M_{\odot})$"

    xmin, xmax, ymin, ymax = 8, 13, 5, 11
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymin + 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'z=0', fontsize=12)

    #Predicted BH-bulge mass relation
    ind = np.where(BH[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = BH[0,0,ind]
    errdn = BH[0,1,ind]
    errup = BH[0,2,ind]
    ax.plot(xplot,yplot[0],color='k',label="Shark (All)")
    ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)

    ind = np.where(BH_morpho[0,0,0,:] != 0)
    xplot = xmf[ind]
    yplot = BH_morpho[0,0,0,ind]
    errdn = BH_morpho[0,0,1,ind]
    errup = BH_morpho[0,0,2,ind]
    ax.plot(xplot,yplot[0],color='Crimson',label="Shark (ETGs)")
    ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='Crimson', alpha=0.25, interpolate=True)

    ind = np.where(BH_morpho[0,1,0,:] != 0)
    xplot = xmf[ind]
    yplot = BH_morpho[0,1,0,ind]
    errdn = BH_morpho[0,1,1,ind]
    errup = BH_morpho[0,1,2,ind]
    ax.plot(xplot,yplot[0],color='Navy',label="Shark (LTGs)")
    ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='Navy', alpha=0.25, interpolate=True)

    #MBH_othermodels = common.load_observation(obsdir, 'Models/SharkVariations/BHBulgeRelation_OtherModels.dat', [0])
    #MBH_f_smbh0p008   = MBH_othermodels[0:29]
    #MBH_f_smbh0p00008 = MBH_othermodels[30:60]
    #ind = np.where(MBH_f_smbh0p008 != 0)
    #xplot = xmf[ind]
    #yplot = MBH_f_smbh0p008[ind]
    #ax.plot(xplot,yplot,color='Goldenrod',linestyle='dashed',label='$f_{\\rm smbh}=8 \\times 10^{-3}$')
    #ind = np.where(MBH_f_smbh0p00008 != 0)
    #xplot = xmf[ind]
    #yplot = MBH_f_smbh0p00008[ind]
    #ax.plot(xplot,yplot,color='Orange',linestyle='dotted',label='$f_{\\rm smbh}=8 \\times 10^{-5}$')

    #BH-bulge relation
    mBH_M13, errup_M13, errdn_M13, mBH_power, mbulge_M13 = common.load_observation(obsdir, 'BHs/MBH_sigma_Mbulge_McConnelMa2013.dat', [0,1,2,3,7])

    ind = np.where((mBH_M13 > 0) & (mbulge_M13 > 0))
    xobs = np.log10(mbulge_M13[ind])
    yobs = np.log10(mBH_M13[ind] * pow(10.0,mBH_power[ind]))
    lerr = np.log10((mBH_M13[ind] - errdn_M13[ind]) * pow(10.0,mBH_power[ind]))
    herr = np.log10((mBH_M13[ind] + errup_M13[ind]) * pow(10.0,mBH_power[ind]))
    ax.errorbar(xobs, yobs, yerr=[yobs-lerr,herr-yobs], ls='None', mfc='None', ecolor = 'r', mec='r',marker='^',label="M13")

    #BH-bulge relation
    #mBH_H04, errup_H04, errdn_H04, mbulge_H04 = common.load_observation(obsdir, 'BHs/MBH_sigma_Mbulge_HaeringRix2004.dat', [0,1,2,4])
    #xobs = np.log10(mbulge_H04)
    #yobs = xobs*0. - 999.
    #indx = np.where( mBH_H04 > 0)
    #yobs[indx] = np.log10(mBH_H04[indx])
    #lerr = yobs*0. - 999.
    #indx = np.where( (mBH_H04-errdn_H04) > 0)
    #lerr[indx]  = np.log10(mBH_H04[indx] - errdn_H04[indx])
    #herr = yobs*0. + 999.
    #indx = np.where( (mBH_H04+errup_H04) > 0)
    #herr[indx]  = np.log10(mBH_H04[indx] + errup_H04[indx])
    #ax.errorbar(xobs, yobs, yerr=[yobs-lerr,herr-yobs], ls='None', mfc='None', ecolor = 'maroon', mec='maroon',marker='s',label="Haering+2004")


    mBH_D19, errup_D19, errdn_D19, mbulge_D19, errbul_D19, mgal_D19, errgal_D19 = common.load_observation(obsdir, 'BHs/MBH_Mstar_Davis2019.dat', [0,1,2,3,4,5,6])
    ind = np.where((mBH_D19 != 0) & (mbulge_D19 != 0))
    ax.errorbar(mbulge_D19[ind], mBH_D19[ind], xerr=errbul_D19[ind], yerr=[errdn_D19[ind], errup_D19[ind]], ls='None', mfc='None', ecolor = 'b', mec='b',marker='o',label="D19 (LTGs)")
    mBH_S19, errbh_S19, mbulge_S19, errbul_S19, mgal_S19, errgal_S19 = common.load_observation(obsdir, 'BHs/MBH_Mstar_Sahu2019.dat', [0,1,2,3,4, 5])
    ind = np.where((mBH_S19 != 0) & (mbulge_S19 != 0))
    ax.errorbar(mbulge_S19[ind], mBH_S19[ind], xerr=errbul_S19[ind], yerr=errbh_S19[ind], ls='None', mfc='None', ecolor = 'Orange', mec='Orange',marker='d',label="S19 (ETGs)")


    common.prepare_legend(ax, ['k','crimson','navy','r','b', 'Orange'], loc=2)
    common.savefig(outdir, fig, 'bulge-BH.pdf')

    #stellar mass-black hole mass relation
    fig = plt.figure(figsize=(6,5.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    xmin, xmax, ymin, ymax = 9, 12.5, 5.5, 11
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymin + 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'z=0', fontsize=12)

    #Predicted BH-stellar mass relation
    ind = np.where(BHSM[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = BHSM[0,0,ind]
    errdn = BHSM[0,1,ind]
    errup = BHSM[0,2,ind]
    ax.plot(xplot,yplot[0],color='k',label="Shark v2.0 (All)")
    ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
    ind = np.where(BHSM_morpho[0,0,0,:] != 0)
    xplot = xmf[ind]
    yplot = BHSM_morpho[0,0,0,ind]
    errdn = BHSM_morpho[0,0,1,ind]
    errup = BHSM_morpho[0,0,2,ind]
    ax.plot(xplot,yplot[0],color='Crimson',label="Shark v2.0 (ETGs)")
    ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='Crimson', alpha=0.25, interpolate=True)

    ind = np.where(BHSM_morpho[0,1,0,:] != 0)
    xplot = xmf[ind]
    yplot = BHSM_morpho[0,1,0,ind]
    errdn = BHSM_morpho[0,1,1,ind]
    errup = BHSM_morpho[0,1,2,ind]
    ax.plot(xplot,yplot[0],color='Navy',label="Shark v2.0 (LTGs)")
    ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='Navy', alpha=0.25, interpolate=True)

    msL18, bhL18_ETGs, bhL18_LTGs = common.load_observation(obsdir, 'Models/SharkVariations/StellarMassBHMass_L18.dat', [0,1,2])
    ind = np.where(bhL18_ETGs > 0)
    ax.plot(msL18[ind], bhL18_ETGs[ind], color='Crimson', linestyle='dashed', label="Shark v1.1 (ETGs)")
    ind = np.where(bhL18_LTGs > 0)
    ax.plot(msL18[ind], bhL18_LTGs[ind], color='Navy', linestyle='dashed', label="Shark v1.1 (LTGs)")


    ind = np.where((mBH_S19 != 0) & (mgal_S19 != 0))
    ax.errorbar(mgal_S19[ind], mBH_S19[ind], xerr=errgal_S19[ind], yerr=errbh_S19[ind], ls='None', mfc='None', ecolor = 'Orange', mec='Orange',marker='d',label="S19 (ETGs)")
    ind = np.where((mBH_D19 != 0) & (mgal_D19 != 0))
    ax.errorbar(mgal_D19[ind], mBH_D19[ind], xerr=errgal_D19[ind], yerr=[errdn_D19[ind],errup_D19[ind]], ls='None', mfc='None', ecolor = 'b', mec='b',marker='d',label="D19 (LTGs)")

    ms, sfr, upperlimflag, mbh, mbherr = common.load_observation(obsdir, 'BHs/MBH_host_gals_Terrazas17.dat', [0,1,2,3,4])
    ind=np.where(sfr-ms > -11.5)
    ax.errorbar(ms[ind], mbh[ind], yerr=mbherr[ind], xerr=0.2, ls='None', mfc='None', ecolor = 'PowderBlue', mec='PowderBlue',marker='s', label="T17 (SF)")
    ind=np.where(sfr-ms <= -11.5)
    ax.errorbar(ms[ind], mbh[ind], yerr=mbherr[ind], xerr=0.2, ls='None', mfc='None', ecolor = 'LightSalmon', mec='LightSalmon',marker='s',label="T17 (P)")

    common.prepare_legend(ax, ['k','crimson','navy', 'crimson','navy', 'Orange', 'b', 'PowderBlue','LightSalmon'], loc=2)
    common.savefig(outdir, fig, 'stellarmass-BH.pdf')

    #SSFR vs BH mass
    fig = plt.figure(figsize=(6,5.5))
    ytit = "$\\rm log_{10} (\\rm sSFR/M_{\odot} yr^{-1})$"
    xtit = "$\\rm log_{10} (\\rm M_{\\rm BH}/M_{\odot})$"

    xmin, xmax, ymin, ymax = 5, 11, -14, -9
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'z=0', fontsize=12)
    ax.text(xleg-0.25, yleg-0.3, '$M_{\\star}> 10^{10}\\, M_{\\odot}$', fontsize=12)

    xL18, yL18, yl_L18, yu_L18 = common.load_observation(obsdir, 'Models/SharkVariations/BHSSFR_Lagos18.dat', [0,1,2,3])
    ax.plot(xL18, yL18,color='k',label="Shark v1.1 (L18)")
    ax.fill_between(xL18,yl_L18,yu_L18, facecolor='k', alpha=0.25, interpolate=True)

    #Predicted BH-bulge mass relation
    ind = np.where(BHSFR[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = BHSFR[0,0,ind]
        errdn = BHSFR[0,1,ind]
        errup = BHSFR[0,2,ind]
        #print("Will print the BH-SSFR correlation")
        #for a,b,c,d in zip(xplot, yplot[0], yplot[0]-errdn[0], yplot[0]+errup[0]):
        #    print(a,b,c,d)
        ax.plot(xplot,yplot[0],color='red',lw=3.5,label="Shark v2.0")
        ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='r', alpha=0.5, interpolate=True)

    #BH-SSFR relation
    ax.errorbar(mbh, sfr-ms, xerr=mbherr, yerr=0.3, ls='None', mfc='None', ecolor = 'b', mec='b',marker='s',label="Terrazas+17")
    ind = np.where(upperlimflag == 1)
    for a,b in zip (mbh[ind], sfr[ind]-ms[ind]):
        ax.arrow(a, b, 0, -0.3, head_width=0.05, head_length=0.1, fc='b', ec='b')

    common.prepare_legend(ax, ['k','r','b'], loc=3)
    common.savefig(outdir, fig, 'BH-SSFR.pdf')

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
        #for i in zip(BT_fractions[0,:]):
        #    print i

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
    ind = np.where(BT_stable1 >= 0)
    xplot = xmf[ind]
    yplot = BT_stable1[ind]
    ax.plot(xplot,yplot,color='Orange',linestyle='dotted',label='$\\epsilon_{\\rm disk}=1$')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    mM16, fM16, errdnfM16, errupfM16 = common.load_observation(obsdir, 'Morph/Moffet16.dat', [0,1,2,3])
    errdnfM16 = np.abs(errdnfM16-fM16)
    errupfM16 = np.abs(errupfM16-fM16)
    ax.errorbar(mM16,fM16,yerr=[errdnfM16,errupfM16], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^',label="Moffett+16")

    common.prepare_legend(ax, ['k', 'r','Goldenrod', 'Orange','grey'], loc=2)
    common.savefig(outdir, fig, 'BTfractions.pdf')

def plot_age_disk(plt, outdir, obsdir, mdisk_z0, rdisk_z0, age_z0):


    fig = plt.figure(figsize=(5,4.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm disk}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm r_{\\star,disk}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.5, 2

    # LTG ##################################
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    im = ax.hexbin(mdisk_z0, rdisk_z0, age_z0, xscale='linear', yscale='linear', gridsize=(20,20), cmap='magma', mincnt=30)
    cbar_ax = fig.add_axes([0.86, 0.15, 0.025, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.set_ylabel('mean stellar age')

    y1 = 0.35 * (xmf - 10) - 0.06 * (1.0) + 1.01
    y4 = 0.35 * (xmf - 10) - 0.06 * (4.0) + 1.01
    y8 = 0.35 * (xmf - 10) - 0.06 * (8.0) + 1.01
    y12 = 0.35 * (xmf - 10) - 0.06 * (12.0) + 1.01

    ax.plot(xmf, y1, linestyle='dotted', color='k')
    ax.plot(xmf, y4, linestyle='dotted', color='k')
    ax.plot(xmf, y8, linestyle='dotted', color='k')
    ax.plot(xmf, y12, linestyle='dotted', color='k')

    common.savefig(outdir, fig, 'age_disks.pdf')

def plot_age_stellar_mass(plt, outdir, obsdir, age_stellar_mass):


    fig = plt.figure(figsize=(5,4.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm t_{\\star}/Gyr$"
    xmin, xmax, ymin, ymax = 8, 12, 0, 14

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    ind = np.where(age_stellar_mass[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = age_stellar_mass[0,0,ind]
        errdn = age_stellar_mass[0,1,ind]
        errup = age_stellar_mass[0,2,ind]
        ax.plot(xplot,yplot[0],color='red',lw=3.5,label="Shark v2.0")
        ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='r', alpha=0.5, interpolate=True)
   

    common.savefig(outdir, fig, 'age_stellar_mass_z0.pdf')


def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'm_bh', 'rstar_disk', 'rstar_bulge', 'type', 
                           'specific_angular_momentum_disk_star', 'specific_angular_momentum_bulge_star',
                           'specific_angular_momentum_disk_gas', 'specific_angular_momentum_bulge_gas',
                           'specific_angular_momentum_disk_gas_atom', 'specific_angular_momentum_disk_gas_mol',
                           'lambda_subhalo', 'mvir_subhalo', 'mgas_disk', 'mgas_bulge','matom_disk', 'mmol_disk', 
                           'matom_bulge', 'mmol_bulge', 'bh_accretion_rate_hh', 'bh_accretion_rate_sb', 'mean_stellar_age', 
                           'sfr_disk', 'sfr_burst')}

    # Loop over redshift and subvolumes
    rcomb = np.zeros(shape = (len(zlist), 3, len(xmf)))
    disk_size = np.zeros(shape = (len(zlist), 3, len(xmf)))
    bulge_size = np.zeros(shape = (len(zlist), 3, len(xmf)))
    bulge_size_mergers = np.zeros(shape = (len(zlist), 3, len(xmf)))
    bulge_size_diskins = np.zeros(shape = (len(zlist), 3, len(xmf)))

    BH = np.zeros(shape = (len(zlist), 3, len(xmf)))
    BHSM = np.zeros(shape = (len(zlist), 3, len(xmf)))
    BH_morpho = np.zeros(shape = (len(zlist), 2, 3, len(xmf)))
    BHSM_morpho = np.zeros(shape = (len(zlist), 2, 3, len(xmf)))

    BHSFR = np.zeros(shape = (len(zlist), 3, len(xmf)))

    disk_size_sat = np.zeros(shape = (len(zlist), 3, len(xmf)))
    disk_size_cen = np.zeros(shape = (len(zlist), 3, len(xmf)))
    BT_fractions = np.zeros(shape = (len(zlist), len(xmf)))
    BT_fractions_nodiskins = np.zeros(shape = (len(zlist), len(xmf)))
    BT_fractions_centrals = np.zeros(shape = (len(zlist), len(xmf)))
    BT_fractions_satellites = np.zeros(shape = (len(zlist), len(xmf)))
    disk_vel =  np.zeros(shape = (len(zlist), 3, len(xmf))) 
    bulge_vel =  np.zeros(shape = (len(zlist), 3, len(xmf)))
    baryonic_TF =  np.zeros(shape = (len(zlist), 3, len(xv))) 
    age_stellar_mass = np.zeros(shape = (len(zlist), 3, len(xmf)))

    bs_error = np.zeros(shape = (len(zlist), len(xmf))) 
    
    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        (mdisk, rdisk, age) = prepare_data(hdf5_data, index, rcomb, disk_size, bulge_size, bulge_size_mergers, bulge_size_diskins, BH,
                     disk_size_sat, disk_size_cen, BT_fractions, BT_fractions_nodiskins, bulge_vel, disk_vel, 
                     BT_fractions_centrals, BT_fractions_satellites, baryonic_TF, BHSM, xmf, xv, bs_error, BHSFR, BH_morpho, BHSM_morpho, 
                     age_stellar_mass)
        if(index == 0):
           mdisk_z0 = mdisk
           rdisk_z0 = rdisk
           age_z0 = age

    plot_sizes(plt, outdir, obsdir, disk_size_cen, disk_size_sat, bulge_size, bulge_size_mergers, bulge_size_diskins)
    plot_velocities(plt, outdir, disk_vel, bulge_vel, baryonic_TF)
    plot_sizes_combined(plt, outdir, obsdir, rcomb)
    plot_bulge_BH(plt, outdir, obsdir, BH, BHSM, BHSFR, BH_morpho, BHSM_morpho)
    plot_bt_fractions(plt, outdir, obsdir, BT_fractions, BT_fractions_nodiskins, BT_fractions_centrals, BT_fractions_satellites)
    plot_age_disk(plt, outdir, obsdir, mdisk_z0, rdisk_z0, age_z0)
    plot_age_stellar_mass(plt, outdir, obsdir, age_stellar_mass)

if __name__ == '__main__':
    main(*common.parse_args())
