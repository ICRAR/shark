#
#    ICRAR - International Centre for Radio Astronomy Research
#    (c) UWA - The University of Western Australia, 2018
#    Copyright by UWA (in the framework of the ICRAR)
#    All rights reserved
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#    MA 02111-1307  USA
#
from scipy.weave.ext_tools import indent
"""Size plots"""

import functools

import numpy as np

import common
import utilities_statistics as us

# Initialize arguments
zlist = ["199","174", "156", "131"]

##################################
#Constants
RExp = 1.67
MpcToKpc = 1e3
G    = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun

mlow = 8
mupp = 12
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0


def prepare_data(hdf5_data, index, rcomb, disk_size, bulge_size, BH,
                 disk_size_sat, disk_size_cen, BT_fractions, bulge_vel, 
                 disk_vel, sam_stars_disk, sam_gas_disk_atom, 
                 sam_gas_disk_mol, sam_halo):

    h0, _, mdisk, mbulge, mBH, rdisk, rbulge, typeg, specific_angular_momentum_disk_star, specific_angular_momentum_bulge_star, specific_angular_momentum_disk_gas, specific_angular_momentum_bulge_gas, specific_angular_momentum_disk_gas_atom, specific_angular_momentum_disk_gas_mol, lambda_sub, mvir_s = hdf5_data
    
    zero_bulge = np.where(rbulge <= 0)
    if(len(rbulge) == len(rbulge[zero_bulge])):
            #case where there is zero bulge build up.
            rbulge[zero_bulge] = 1.0
            specific_angular_momentum_bulge_star[zero_bulge] = 1.0
            mbulge[zero_bulge] = 10.0

    bin_it = functools.partial(us.wmedians, xbins=xmf)

    sam_subhalo = 1.41421356237 * G**0.66 * lambda_sub * mvir_s**0.66 / (h0*100.0)**0.33;

    vdisk = specific_angular_momentum_disk_star / rdisk / 2.0 #in km/s
    vbulge = specific_angular_momentum_bulge_star / rbulge / 2.0 #in km/s
    
    ind = np.where(mdisk+mbulge > 0)
    rcomb[index,:] = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                            y=np.log10((mdisk[ind]*rdisk[ind]  + mbulge[ind]*rbulge[ind])*MpcToKpc / (mdisk[ind]+mbulge[ind])))
    BT_fractions[index] = us.fractional_contribution(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),y=mbulge[ind]/(mdisk[ind]+mbulge[ind]), xbins=xmf)

    ind = np.where(mdisk > 0)
    disk_size[index,:] = bin_it(x=np.log10(mdisk[ind]) - np.log10(float(h0)),
                                y=np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)))

    ind = np.where((specific_angular_momentum_disk_star > 0) & (mdisk+mbulge > 0) & (mbulge/mdisk < 0.5) & (typeg == 0))
    sam_stars_disk[index,:] = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                y=np.log10(specific_angular_momentum_disk_star[ind]) - np.log10(float(h0)))
    
    ind = np.where((specific_angular_momentum_disk_gas > 0) & (mdisk+mbulge > 0) & (mbulge/mdisk < 0.5) & (typeg == 0))
    sam_gas_disk_atom[index,:]   = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                y=np.log10(specific_angular_momentum_disk_gas[ind]) - np.log10(float(h0)))
    
    ind = np.where((specific_angular_momentum_disk_gas_mol > 0) & (mdisk+mbulge > 0) & (mbulge/mdisk < 0.5) & (typeg == 0))
    sam_gas_disk_mol[index,:]   = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                y=np.log10(specific_angular_momentum_disk_gas_mol[ind]) - np.log10(float(h0)))
    
    ind = np.where((sam_subhalo > 0) & (mdisk+mbulge > 0) & (mbulge/mdisk < 0.5) & (typeg == 0))
    sam_halo[index,:]       = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                y=np.log10(sam_subhalo[ind]) - np.log10(float(h0)))

    ind = np.where((mdisk > 0) & (typeg == 0) & (mbulge/mdisk < 0.5))
    disk_vel[index,:] = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                y=np.log10(vdisk[ind]))

    ind = np.where((mdisk > 0) & (typeg == 0))
    disk_size_cen[index,:]  = bin_it(x=np.log10(mdisk[ind]) - np.log10(float(h0)),
                                    y=np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)))

    ind = np.where((mdisk > 0) & (typeg > 0))
    disk_size_sat[index,:] = bin_it(x=np.log10(mdisk[ind]) - np.log10(float(h0)),
                                    y=np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)))

    ind = np.where(mbulge > 0)
    bulge_size[index,:] = bin_it(x=np.log10(mbulge[ind]) - np.log10(float(h0)),
                                 y=np.log10(rbulge[ind]*MpcToKpc) - np.log10(float(h0)))

    BH[index,:] = bin_it(x=np.log10(mbulge[ind]) - np.log10(float(h0)),
                    y=np.log10(mBH[ind]) - np.log10(float(h0)))
    
    ind = np.where((mbulge > 0) & (mbulge/mdisk > 0.5))
    bulge_vel[index,:] = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                    y=np.log10(vbulge[ind]))
    


def plot_sizes(plt, outdir, obsdir, disk_size_cen, disk_size_sat, bulge_size):

    fig = plt.figure(figsize=(5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm stars,disk}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm r_{\\rm 50, disk}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.1, 2
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # LTG ##################################
    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'z=0')

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(disk_size_cen[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = disk_size_cen[0,0,ind]
    errdn = disk_size_cen[0,1,ind]
    errup = disk_size_cen[0,2,ind]
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="SHArk centrals")

    #Predicted size-mass for disks in disk=dominated galaxies satellites
    ind = np.where(disk_size_sat[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = disk_size_sat[0,0,ind]
    errdn = disk_size_sat[0,1,ind]
    errup = disk_size_sat[0,2,ind]
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'r', mec='r',marker='v',markersize='5',label="SHArk satellites")

    #Lange et al. (2016)
    a = 5.56
    aerr = 1.745
    b = 0.274
    ind = np.where(xmf < 10.3)
    rL16 = np.log10(a*pow((pow(10.0,xmf)/1e10),b))

    ax.plot(xmf[ind],rL16[ind],'b', linestyle='solid', label ='L16 disks')


    m,r = common.load_observation(obsdir, 'rdisk_L16.dat', [0,1])
    ax.plot(m[0:36], r[0:36], linestyle='dotted',color='b',label="50th, 68th, 90th")
    ax.plot(m[38:83], r[38:83], linestyle='dotted',color='b')
    ax.plot(m[85:128], r[85:129], linestyle='dotted',color='b')

    common.prepare_legend(ax, ['b','b','k','r'], loc=2)

    # ETGs ##################################
    xtit = "$\\rm log_{10} (\\rm M_{\\rm stars, bulge}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm r_{\\rm 50, bulge}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.5, 2
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(212)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'z=0')

    #Predicted size-mass for bulges in bulge-dominated systems
    ind = np.where(bulge_size[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = bulge_size[0,0,ind]
        errdn = bulge_size[0,1,ind]
        errup = bulge_size[0,2,ind]
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="SHArk bulges")

    #Lange et al. (2016)
    a = 2.319
    aerr = 1.186
    b = 0.19565217391
    a2 = 1.390
    aerr2 = 0.9
    b2 = 0.624

    ind = np.where(xmf <= 10.3)

    rL16 = np.log10(a*pow((pow(10.0,xmf[ind])/1e10),b))
    ax.plot(xmf[ind],rL16,'r', linestyle='solid', label ='L16 $M_{\\rm stars}<2\\times 10^{10}\\rm M_{\odot}$')

    ind = np.where((xmf >= 10.3) & (xmf < 11.5))
    rL16_2 = np.log10(a2*pow((pow(10.0,xmf[ind])/1e10),b2))

    ax.plot(xmf[ind],rL16_2,'m', linestyle='solid', label ='L16 $M_{\\rm stars}>2\\times 10^{10}\\rm M_{\odot}$')
    
    m,r = common.load_observation(obsdir, 'rbulge_L16.dat', [0,1])
    ax.plot(m[0:39], r[0:39], linestyle='dotted',color='r',label="50th, 68th, 90th")
    ax.plot(m[41:76], r[41:76], linestyle='dotted',color='r')
    ax.plot(m[78:115], r[78:115], linestyle='dotted',color='r')

    common.prepare_legend(ax, ['r','m','r','k'], loc=2)
    common.savefig(outdir, fig, 'sizes.pdf')


def plot_velocities(plt, outdir, disk_vel, bulge_vel):

    fig = plt.figure(figsize=(5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm stars}/M_{\odot})$"
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
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="SHArk centrals B/T<0.5")

    common.prepare_legend(ax, ['k'], loc=2)

    # ETGs ##################################
    xtit = "$\\rm log_{10} (\\rm M_{\\rm stars}/M_{\odot})$"
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
        ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="SHArk bulges B/T > 0.5")

        common.prepare_legend(ax, ['k'], loc=2)

    common.savefig(outdir, fig, 'velocities.pdf')
 
def plot_specific_am(plt, outdir, obsdir, sam_stars_disk, sam_gas_disk_atom, sam_gas_disk_mol, sam_halo):

    fig = plt.figure(figsize=(5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm stars}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm j_{\\rm disk}/kpc\\, km s^{-1})$"
    xmin, xmax, ymin, ymax = 8, 12, 1.5, 5
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # LTG ##################################
    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'z=0')

    #Romanowsky and Fall (2012)
    
    #a = 0.67
    #b = -4.029
    #berr = 0.2
    #ind = np.where(xmf < 11)
    #rRF12 = a * xmf[ind] + b
    #rRF12Low  = a * xmf[ind] + b - berr
    #rRF12High = a * xmf[ind] + b + berr

    #ax.plot(xmf[ind],rRF12,'g', linestyle='solid', label ='R&F')
    #ax.plot(xmf[ind],rRF12Low,'g', linestyle='dotted')
    #ax.plot(xmf[ind],rRF12High,'g', linestyle='dotted')

    ind = np.where(sam_halo[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = sam_halo[0,0,ind] + 3.0
    errdn = sam_halo[0,1,ind]
    errup = sam_halo[0,2,ind]
    ax.plot(xplot,yplot[0],color='k',label="DM")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.2,interpolate=True)

    #Predicted sAM-mass for disks in disk=dominated galaxies
    ind = np.where(sam_stars_disk[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = sam_stars_disk[0,0,ind]+ 3.0
    errdn = sam_stars_disk[0,1,ind]
    errup = sam_stars_disk[0,2,ind]
    ax.plot(xplot,yplot[0],color='r',label="stars")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.5,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.5,interpolate=True)

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(sam_gas_disk_atom[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = sam_gas_disk_atom[0,0,ind] + 3.0
    errdn = sam_gas_disk_atom[0,1,ind]
    errup = sam_gas_disk_atom[0,2,ind]
    ax.plot(xplot,yplot[0],color='b',label="atomic ISM")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.5,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.5,interpolate=True)
    
    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(sam_gas_disk_mol[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = sam_gas_disk_mol[0,0,ind] + 3.0
    errdn = sam_gas_disk_mol[0,1,ind]
    errup = sam_gas_disk_mol[0,2,ind]
    ax.plot(xplot,yplot[0],color='g',label="molecular ISM")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', alpha=0.5,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', alpha=0.5,interpolate=True)

    bt, ms, mg, js, jg, jmol = common.load_observation(obsdir, 'Obreschkow14_FP.dat', [2,7,8,12,14,15])
    ind = np.where(bt < 0.5)
    ax.plot(ms[ind], js[ind], 'ro',fillstyle='full', label="Obreschkow+14")
    ax.plot(ms[ind], jg[ind], 'bo',fillstyle='full')
    ax.plot(ms[ind], jmol[ind], 'go',fillstyle='full')
    
    mg, ms, jg, js = common.load_observation(obsdir, 'LITTLETHINGS_Butler16.dat', [1,3,7,9])
    ax.plot(ms, js, 'rs',fillstyle='full',label="Butler+16")
    ax.plot(ms, jg, 'bs',fillstyle='full')

    common.prepare_legend(ax, ['k'], loc=2)

    ###################################
    ax = fig.add_subplot(212)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'z=2')

    ind = np.where(sam_halo[3,0,:] != 0)
    xplot = xmf[ind]
    yplot = sam_halo[3,0,ind] + 3.0
    errdn = sam_halo[3,1,ind]
    errup = sam_halo[3,2,ind]
    ax.plot(xplot,yplot[0],color='k',label="DM")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.2,interpolate=True)

    #Predicted sAM-mass for disks in disk=dominated galaxies
    ind = np.where(sam_stars_disk[3,0,:] != 0)
    xplot = xmf[ind]
    yplot = sam_stars_disk[3,0,ind]+ 3.0
    errdn = sam_stars_disk[3,1,ind]
    errup = sam_stars_disk[3,2,ind]
    ax.plot(xplot,yplot[0],color='r',label="stars")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.5,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.5,interpolate=True)

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(sam_gas_disk_atom[3,0,:] != 0)
    xplot = xmf[ind]
    yplot = sam_gas_disk_atom[3,0,ind] + 3.0
    errdn = sam_gas_disk_atom[3,1,ind]
    errup = sam_gas_disk_atom[3,2,ind]
    ax.plot(xplot,yplot[0],color='b',label="atomic ISM")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.5,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.5,interpolate=True)

    ind = np.where(sam_gas_disk_mol[3,0,:] != 0)
    xplot = xmf[ind]
    yplot = sam_gas_disk_mol[3,0,ind] + 3.0
    errdn = sam_gas_disk_mol[3,1,ind]
    errup = sam_gas_disk_mol[3,2,ind]
    ax.plot(xplot,yplot[0],color='g',label="molecular ISM")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', alpha=0.5,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', alpha=0.5,interpolate=True)

    common.prepare_legend(ax, ['k'], loc=2)

    common.savefig(outdir, fig, 'specific_am.pdf')
    
def plot_sizes_combined(plt, outdir, rcomb):

    fig = plt.figure(figsize=(5,5))

    # Total ##################################
    xtit="$\\rm log_{10} (\\rm M_{\\rm stars, total}/M_{\odot})$"
    ytit="$\\rm log_{10} (\\rm r_{\\rm 50, comb}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.5, 2

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    #Predicted size-mass for disks
    ind = np.where(rcomb[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = rcomb[0,0,ind]
    errdn = rcomb[0,1,ind]
    errup = rcomb[0,2,ind]
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="SHArk disk+bulge combined")

    common.prepare_legend(ax, ['k','k','k'], loc=2)
    common.savefig(outdir, fig, 'sizes_combined.pdf')


def plot_bulge_BH(plt, outdir, obsdir, BH):

    fig = plt.figure(figsize=(5,5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm bulge}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm M_{\\rm BH}/M_{\odot})$"

    xmin, xmax, ymin, ymax = 8, 12, 5, 11
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'z=0')

    #Predicted SMHM
    ind = np.where(BH[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = BH[0,0,ind]
        errdn = BH[0,1,ind]
        errup = BH[0,2,ind]
        ax.plot(xplot,yplot[0],color='k',label="SHArk")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='grey', interpolate=True)

    #BH-bulge relation
    mBH_H04, errup_H04, errdn_H04, mbulge_H04 = common.load_observation(obsdir, 'MBH_sigma_Mbulge_HaeringRix2004.data', [0,1,2,4])

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
    ax.errorbar(xobs, yobs, yerr=[yobs-lerr,herr-yobs], ls='None', mfc='None', ecolor = 'r', mec='r',marker='^',label="Haering+04")

    common.prepare_legend(ax, ['k','r'], loc=2)
    common.savefig(outdir, fig, 'bulge-BH.pdf')


def plot_bt_fractions(plt, outdir, obsdir, BT_fractions):

    fig = plt.figure(figsize=(5,5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm stars}/M_{\odot})$"
    ytit = "$\\rm f_{\\rm bulges}$"
    xmin, xmax, ymin, ymax = 8, 12, -0.1, 1.05

    # LTG ##################################
    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #Predicted size-mass for disks
    ind = np.where(BT_fractions[0,:] >= 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = BT_fractions[0,ind]
        ax.plot(xplot,yplot[0],'r', label ='all z=0 galaxies')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    mM16, fM16, errdnfM16, errupfM16 = common.load_observation(obsdir, 'Moffet16.dat', [0,1,2,3])
    errdnfM16 = np.abs(errdnfM16-fM16)
    errupfM16 = np.abs(errupfM16-fM16)
    ax.errorbar(mM16,fM16,yerr=[errdnfM16,errupfM16], ls='None', mfc='None', ecolor = 'b', mec='b',marker='^',label="Moffett+16")

    common.prepare_legend(ax, ['r', 'b'], loc=2)
    common.savefig(outdir, fig, 'BTfractions.pdf')


def main():

    plt = common.load_matplotlib()
    fields = {'Galaxies': ('mstars_disk', 'mstars_bulge', 'mBH',
                           'rdisk_star', 'rbulge_star', 'type', 
                           'specific_angular_momentum_disk_star', 'specific_angular_momentum_bulge_star',
                           'specific_angular_momentum_disk_gas', 'specific_angular_momentum_bulge_gas',
                           'specific_angular_momentum_disk_gas_atom', 'specific_angular_momentum_disk_gas_mol',
                           'lambda_subhalo', 'mvir_subhalo')}

    modeldir, outdir, obsdir = common.parse_args(requires_snapshot=False)

    # Loop over redshift and subvolumes
    rcomb = np.zeros(shape = (len(zlist), 3, len(xmf)))
    disk_size = np.zeros(shape = (len(zlist), 3, len(xmf)))
    bulge_size = np.zeros(shape = (len(zlist), 3, len(xmf)))
    BH = np.zeros(shape = (len(zlist), 3, len(xmf)))
    disk_size_sat = np.zeros(shape = (len(zlist), 3, len(xmf)))
    disk_size_cen = np.zeros(shape = (len(zlist), 3, len(xmf)))
    BT_fractions = np.zeros(shape = (len(zlist), len(xmf)))
    disk_vel =  np.zeros(shape = (len(zlist), 3, len(xmf))) 
    bulge_vel =  np.zeros(shape = (len(zlist), 3, len(xmf)))
    sam_stars_disk    = np.zeros(shape = (len(zlist), 3, len(xmf)))
    sam_gas_disk_atom = np.zeros(shape = (len(zlist), 3, len(xmf)))
    sam_gas_disk_mol  = np.zeros(shape = (len(zlist), 3, len(xmf)))
    
    sam_halo       = np.zeros(shape = (len(zlist), 3, len(xmf)))

    for index in range(0,4):
        hdf5_data = common.read_data(modeldir, zlist[index], fields, subvolume=0)
        prepare_data(hdf5_data, index, rcomb, disk_size, bulge_size, BH,
                     disk_size_sat, disk_size_cen, BT_fractions, bulge_vel, disk_vel, 
                     sam_stars_disk, sam_gas_disk_atom, sam_gas_disk_mol, sam_halo)

    plot_sizes(plt, outdir, obsdir, disk_size_cen, disk_size_sat, bulge_size)
    plot_velocities(plt, outdir, disk_vel, bulge_vel)
    plot_specific_am(plt, outdir, obsdir, sam_stars_disk, sam_gas_disk_atom, sam_gas_disk_mol, sam_halo)
    plot_sizes_combined(plt, outdir, rcomb)
    plot_bulge_BH(plt, outdir, obsdir, BH)
    plot_bt_fractions(plt, outdir, obsdir, BT_fractions)

if __name__ == '__main__':
    main()
