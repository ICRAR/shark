
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
"""Angular momentum plots"""

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
Omegab   = 0.0491
OmegaM   = 0.3121
fbar     = Omegab/(OmegaM-Omegab)

#halo mass bins
mlow = 10.0
mupp = 15.5
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

dmobs = 0.4
mbins_obs = np.arange(mlow,mupp,dmobs)
xmf_obs = mbins_obs + dmobs/2.0

llow = -3.0
lupp = 6.0
dl   = 0.15
lbins = np.arange(llow,lupp,dl)
xlf   = lbins + dl/2.0


def prepare_data(hdf5_data, index, sam_vs_sam_halo_disk, sam_vs_sam_halo_gal,
                 sam_vs_sam_halo_disk_gas, sam_vs_sam_halo_bar, m_vs_m_halo_disk, m_vs_m_halo_gal, 
                 m_vs_m_halo_disk_gas, m_vs_m_halo_bar, r_vs_r_halo_disk, r_vs_r_halo_gal, 
                 r_vs_r_halo_disk_gas, r_vs_r_halo_bar):

    (h0, _, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     mBH, rdisk, rbulge, typeg, specific_angular_momentum_disk_star, specific_angular_momentum_bulge_star, 
     specific_angular_momentum_disk_gas, specific_angular_momentum_bulge_gas, specific_angular_momentum_disk_gas_atom, 
     specific_angular_momentum_disk_gas_mol, lambda_sub, mvir_s, mvir, matom_disk, mmol_disk, mgas_disk,
     matom_bulge, mmol_bulge, mgas_bulge, sfr_disk, sfr_bulge, vvir, rgas_disk, rgas_bulge) = hdf5_data

    rvir  = G * mvir / pow(vvir,2.0) * MpcToKpc / h0
    ind = np.where(rvir == max(rvir)) 
    print 'maximum rvir', np.log10(max(rvir)), mvir[ind]
    ind = np.where(mvir == max(mvir)) 
    print 'maximum mvir', np.log10(rvir[ind]), max(mvir)

    rbar  = (rdisk * mdisk + rbulge * mbulge + rgas_disk * mgas_disk + rgas_bulge * mgas_bulge) / (mdisk + mbulge + mgas_disk + mgas_bulge) * MpcToKpc / h0
    rgal  = (rdisk * mdisk + rbulge * mbulge)/(mdisk + mbulge) * MpcToKpc / h0
    rstar_disk = rdisk * MpcToKpc / h0
    rgas_disk  = rgas_disk * MpcToKpc / h0
 
    mbulge_mergers = mburst_mergers + mstars_bulge_mergers_assembly
    zero_bulge = np.where(rbulge <= 0)
    if(len(rbulge) == len(rbulge[zero_bulge])):
            #case where there is zero bulge build up.
            rbulge[zero_bulge] = 1e-10
            specific_angular_momentum_bulge_star[zero_bulge] = 1.0
            mbulge[zero_bulge] = 10.0

    bin_it = functools.partial(us.wmedians, xbins=xmf, low_numbers=True)
    bin_it_j    = functools.partial(us.wmedians, xbins=xlf, low_numbers=True)

    sam_subhalo = 1.41421356237 * G**0.66 * lambda_sub * mvir_s**0.66 / (h0*100.0)**0.33
    sam_hhalo   = 1.41421356237 * G**0.66 * lambda_sub * mvir**0.66 / (h0*100.0)**0.33

    #calculate effective specific angular momentum by mass-weighting the contributions from the disk and bulge
    jstars      = (specific_angular_momentum_disk_star * mdisk + specific_angular_momentum_bulge_star * mbulge) / (mdisk+mbulge)
    jbar        = (specific_angular_momentum_disk_star * mdisk + specific_angular_momentum_disk_gas * mgas_disk + 
                   specific_angular_momentum_bulge_star * mbulge + specific_angular_momentum_bulge_gas * mgas_bulge) / (mdisk + mbulge + mgas_disk + mgas_bulge)
    mbar        = (mdisk + mbulge + mgas_disk + mgas_bulge)

    vdisk  = specific_angular_momentum_disk_star / rdisk / 2.0 #in km/s
    vbulge = specific_angular_momentum_bulge_star / rbulge / 2.0 #in km/s

    specific_angular_momentum_disk = (specific_angular_momentum_disk_star * mdisk + specific_angular_momentum_disk_gas * mgas_disk) / (mdisk + mgas_disk)

    vr_halo = G**0.66 * mvir**0.66 / (h0*100.0)**0.33
    lh = lambda_sub
    lj = np.zeros(shape = (4, len(mdisk))) 
    lm = np.zeros(shape = (4, len(mdisk))) 
    
    lj[0,:] = specific_angular_momentum_disk  / 1.41421356237 / vr_halo
    lj[1,:] = specific_angular_momentum_disk_star  / 1.41421356237 / vr_halo
    lj[2,:] = specific_angular_momentum_disk_gas   / 1.41421356237 / vr_halo
    lj[3,:] = jstars / 1.41421356237 / vr_halo

    lm[0,:] = specific_angular_momentum_disk  / 1.41421356237 / G**0.66 * (h0*100.0)**0.33 / (mdisk + mgas_disk)**0.66
    lm[1,:] = specific_angular_momentum_disk_star  / 1.41421356237 / G**0.66 * (h0*100.0)**0.33 / (mdisk)**0.66
    lm[2,:] = specific_angular_momentum_disk_gas   / 1.41421356237 / G**0.66 * (h0*100.0)**0.33 / (mgas_disk)**0.66
    lm[3,:] = jstars  / 1.41421356237 / G**0.66 * (h0*100.0)**0.33 / (mdisk + mbulge)**0.66

    bt = np.zeros(shape = (len(mdisk)))
    ms = np.zeros(shape = (len(mdisk)))
    ssfr = np.zeros(shape = (len(mdisk)))
    ind = np.where(mdisk+mbulge > 0) 
    bt[ind] = mbulge[ind] / (mdisk[ind] + mbulge[ind])
    ms[ind] = np.log10(mdisk[ind] + mbulge[ind])
    ind = np.where((mdisk+mbulge > 0) & (sfr_disk + sfr_bulge > 0))
    ssfr[ind] = np.log10((sfr_disk[ind] + sfr_bulge[ind])) - ms[ind] #in Gyr
 
    thresh = [0, 0.5]
    mass_cut = 7.5
    for c in range(0,2):
        ind = np.where((jbar > 0) & (mbar > mass_cut)  & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c]))
        sam_vs_sam_halo_bar[index,:,:,c]      = bin_it_j(x=np.log10(sam_hhalo[ind])  - np.log10(float(h0)),
                                                         y=np.log10(jbar[ind])  - np.log10(float(h0)))
        m_vs_m_halo_bar[index,:,:,c]          = bin_it(x=np.log10(mvir[ind])  - np.log10(float(h0)),
                                                       y=np.log10(mbar[ind])  - np.log10(float(h0)))
        r_vs_r_halo_bar[index,:,:,c]          = bin_it_j(x=np.log10(rvir[ind])  - np.log10(float(h0)),
                                                       y=np.log10(rbar[ind])  - np.log10(float(h0)))

        ind = np.where((jstars > 0) & (ms > mass_cut)  & (mdisk+mbulge > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c]))
        sam_vs_sam_halo_gal[index,:,:,c]      = bin_it_j(x=np.log10(sam_hhalo[ind])  - np.log10(float(h0)),
                                                  y=np.log10(jstars[ind])  - np.log10(float(h0)))

        m_vs_m_halo_gal[index,:,:,c]          = bin_it(x=np.log10(mvir[ind])  - np.log10(float(h0)),
                                                  y=ms[ind]  - np.log10(float(h0)))
        r_vs_r_halo_gal[index,:,:,c]          = bin_it_j(x=np.log10(rvir[ind])  - np.log10(float(h0)),
                                                       y=np.log10(rgal[ind])  - np.log10(float(h0)))

        ind = np.where((specific_angular_momentum_disk_star > 0) & (ms > mass_cut)  & (mdisk+mbulge > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c]))
        sam_vs_sam_halo_disk[index,:,:,c]     = bin_it_j(x=np.log10(sam_hhalo[ind])  - np.log10(float(h0)),
                                                         y=np.log10(specific_angular_momentum_disk_star[ind])  - np.log10(float(h0)))
        m_vs_m_halo_disk[index,:,:,c]         = bin_it(x=np.log10(mvir[ind])  - np.log10(float(h0)),
                                                       y=np.log10(mdisk[ind])  - np.log10(float(h0)))
        r_vs_r_halo_disk[index,:,:,c]         = bin_it_j(x=np.log10(rvir[ind])  - np.log10(float(h0)),
                                                       y=np.log10(rstar_disk[ind])  - np.log10(float(h0)))
  
        ind = np.where((specific_angular_momentum_disk_gas > 0) & (ms > mass_cut)  & (mdisk+mbulge > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c]))
        sam_vs_sam_halo_disk_gas[index,:,:,c] = bin_it_j(x=np.log10(sam_hhalo[ind])  - np.log10(float(h0)),
                                                  y=np.log10(specific_angular_momentum_disk_gas[ind])  - np.log10(float(h0)))
        m_vs_m_halo_disk_gas[index,:,:,c]     = bin_it(x=np.log10(mvir[ind])  - np.log10(float(h0)),
                                                   y=np.log10(mgas_disk[ind])  - np.log10(float(h0)))
        r_vs_r_halo_disk_gas[index,:,:,c]     = bin_it_j(x=np.log10(rvir[ind])  - np.log10(float(h0)),
                                                       y=np.log10(rgas_disk[ind])  - np.log10(float(h0)))
   

    return (lh, lj, lm, bt, ms, ssfr)

def plot_specific_am_ratio(plt, outdir, obsdir, sam_vs_sam_halo_disk, sam_vs_sam_halo_gal, sam_vs_sam_halo_disk_gas, sam_vs_sam_halo_bar, 
        m_vs_m_halo_disk, m_vs_m_halo_gal, m_vs_m_halo_disk_gas, m_vs_m_halo_bar, r_vs_r_halo_disk, 
        r_vs_r_halo_gal, r_vs_r_halo_disk_gas, r_vs_r_halo_bar):

    selec = 0 #all galaxies

    #plot specific AM vs. specific AM 
    fig = plt.figure(figsize=(9,8))
    xtit = "$\\rm log_{10} (\\rm j_{\\rm halo}/kpc\,km\,s^{-1})$"
    ytit = "$\\rm log_{10} (\\rm j_{\\rm gas}/kpc\,km\,s^{-1}$)"
    xmin, xmax, ymin, ymax = 1,6,1,6
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (221, 222, 223, 224)
    indz = (0, 1, 2, 3)
    zinplot = (0, 0.5, 1, 2) 

    # LTG ##################################
    for z,s,p in zip(zinplot, indz, subplots):
	    ax = fig.add_subplot(p)
	    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
            ax.text(xleg, yleg, 'z=%s' % str(z), fontsize=10)
            #if(s == 0):
            #    ax.text(5.5,6.4,'Lagos+18',fontsize=14)

	    ind = np.where(sam_vs_sam_halo_gal[s,0,:,selec] != 0)
	    xplot = xlf[ind]+3
	    yplot = sam_vs_sam_halo_gal[s,0,ind,selec]+3
	    errdn = sam_vs_sam_halo_gal[s,1,ind,selec]
	    errup = sam_vs_sam_halo_gal[s,2,ind,selec]
	    ax.plot(xplot,yplot[0],color='k',label="all stars")
	    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.2,interpolate=True)
	    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.2,interpolate=True)

	    ind = np.where(sam_vs_sam_halo_disk[s,0,:,selec] != 0)
	    xplot = xlf[ind]+3
	    yplot = sam_vs_sam_halo_disk[s,0,ind,selec]+3
	    errdn = sam_vs_sam_halo_disk[s,1,ind,selec]
	    errup = sam_vs_sam_halo_disk[s,2,ind,selec]
	    ax.plot(xplot,yplot[0],color='g',label="disk stars")
	    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', alpha=0.2,interpolate=True)
	    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', alpha=0.2,interpolate=True)

	    ind = np.where(sam_vs_sam_halo_disk_gas[s,0,:,selec] != 0)
	    xplot = xlf[ind]+3
	    yplot = sam_vs_sam_halo_disk_gas[s,0,ind,selec]+3
	    errdn = sam_vs_sam_halo_disk_gas[s,1,ind,selec]
	    errup = sam_vs_sam_halo_disk_gas[s,2,ind,selec]
	    ax.plot(xplot,yplot[0],color='b',label="disk gas")
	    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.2,interpolate=True)
	    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.2,interpolate=True)

            xplot = [1,5]
            ax.plot(xplot,xplot,color='grey',linestyle='dotted')
            if(s == 0):
               common.prepare_legend(ax, ['k','g','b'], loc=2)

    common.savefig(outdir, fig, 'specific_am_halo_vs_galaxy.pdf')

    #plot specific AM vs. specific AM z=0 only
    fig = plt.figure(figsize=(4,11.5))
    xtit = "$\\rm log_{10} (\\rm j_{\\rm halo}/kpc\,km\,s^{-1})$"
    ytit = "$\\rm log_{10} (\\rm j_{\\rm gal}/kpc\,km\,s^{-1}$)"
    xmin, xmax, ymin, ymax = 1,6,1,5
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # LTG ##################################
    ax = fig.add_subplot(311)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    plt.subplots_adjust(left=0.2)
    #ax.text(2.7,5.1,'Shark-AM',fontsize=14)
    #ax.text(2.5,5.1,'Shark-default',fontsize=14)
    #ax.text(0.99,5.1,'No AGN feedback, no disk-inst',fontsize=14)
    ax.text(2.45,5.1,'No disk-inst',fontsize=14)

    s = 0 #z=0
    ind = np.where(sam_vs_sam_halo_bar[s,0,:,selec] != 0)
    xplot = xlf[ind]+3
    yplot = sam_vs_sam_halo_bar[s,0,ind,selec]+3
    errdn = sam_vs_sam_halo_bar[s,1,ind,selec]
    errup = sam_vs_sam_halo_bar[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='k',label="all baryons")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.2,interpolate=True)

    ind = np.where(sam_vs_sam_halo_gal[s,0,:,selec] != 0)
    xplot = xlf[ind]+3
    yplot = sam_vs_sam_halo_gal[s,0,ind,selec]+3
    errdn = sam_vs_sam_halo_gal[s,1,ind,selec]
    errup = sam_vs_sam_halo_gal[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='r',linestyle='dashdot', label="all stars")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.2,interpolate=True)

    ind = np.where(sam_vs_sam_halo_disk[s,0,:,selec] != 0)
    xplot = xlf[ind]+3
    yplot = sam_vs_sam_halo_disk[s,0,ind,selec]+3
    errdn = sam_vs_sam_halo_disk[s,1,ind,selec]
    errup = sam_vs_sam_halo_disk[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='g',linestyle='dashed',label="disk stars")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', alpha=0.2,interpolate=True)

    ind = np.where(sam_vs_sam_halo_disk_gas[s,0,:,selec] != 0)
    xplot = xlf[ind]+3
    yplot = sam_vs_sam_halo_disk_gas[s,0,ind,selec]+3
    errdn = sam_vs_sam_halo_disk_gas[s,1,ind,selec]
    errup = sam_vs_sam_halo_disk_gas[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='b', linestyle='dotted', label="disk gas")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.2,interpolate=True)

    ax.plot(xlf,xlf,color='grey',linestyle='dotted')

    common.prepare_legend(ax, ['k','r','g','b'], loc=2)

    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm M_{\\rm gal}(\\Omega_{\\rm DM}/\\Omega_{\\rm b})/M_{\odot})$"
    xmin, xmax, ymin, ymax = 10,15,8,13.5
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ##################################
    ax = fig.add_subplot(312)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    plt.subplots_adjust(left=0.2)

    ind = np.where(m_vs_m_halo_bar[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = m_vs_m_halo_bar[s,0,ind,selec]-np.log10(fbar)
    print 'fbar',np.log10(fbar),fbar,yplot[0]
    errdn = m_vs_m_halo_bar[s,1,ind,selec]
    errup = m_vs_m_halo_bar[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='k')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.2,interpolate=True)

    ind = np.where(m_vs_m_halo_gal[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = m_vs_m_halo_gal[s,0,ind,selec]-np.log10(fbar)
    errdn = m_vs_m_halo_gal[s,1,ind,selec]
    errup = m_vs_m_halo_gal[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='r',linestyle='dashdot')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.2,interpolate=True)

    ind = np.where(m_vs_m_halo_disk[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = m_vs_m_halo_disk[s,0,ind,selec]-np.log10(fbar)
    errdn = m_vs_m_halo_disk[s,1,ind,selec]
    errup = m_vs_m_halo_disk[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='g',linestyle='dashed')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', alpha=0.2,interpolate=True)

    ind = np.where(m_vs_m_halo_disk_gas[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = m_vs_m_halo_disk_gas[s,0,ind,selec]-np.log10(fbar)
    errdn = m_vs_m_halo_disk_gas[s,1,ind,selec]
    errup = m_vs_m_halo_disk_gas[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='b', linestyle='dotted')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.2,interpolate=True)

    ax.plot(xmf,xmf,color='grey',linestyle='dotted')

    xtit = "$\\rm log_{10} (\\rm r_{\\rm halo}/kpc)$"
    ytit = "$\\rm log_{10} (\\rm r_{\\rm gal}/kpc$)"
    xmin, xmax, ymin, ymax = 1.7,3.3,0,1.5
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ##################################
    ax = fig.add_subplot(313)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    plt.subplots_adjust(left=0.2)
    s = 0 #z=0
    ind = np.where(r_vs_r_halo_bar[s,0,:,selec] != 0)
    xplot = xlf[ind]
    yplot = r_vs_r_halo_bar[s,0,ind,selec]
    errdn = r_vs_r_halo_bar[s,1,ind,selec]
    errup = r_vs_r_halo_bar[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='k')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.2,interpolate=True)

    ind = np.where(r_vs_r_halo_gal[s,0,:,selec] != 0)
    xplot = xlf[ind]
    yplot = r_vs_r_halo_gal[s,0,ind,selec]
    errdn = r_vs_r_halo_gal[s,1,ind,selec]
    errup = r_vs_r_halo_gal[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='r',linestyle='dashdot')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.2,interpolate=True)

    ind = np.where(r_vs_r_halo_disk[s,0,:,selec] != 0)
    xplot = xlf[ind]
    yplot = r_vs_r_halo_disk[s,0,ind,selec]
    errdn = r_vs_r_halo_disk[s,1,ind,selec]
    errup = r_vs_r_halo_disk[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='g',linestyle='dashed')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', alpha=0.2,interpolate=True)

    ind = np.where(r_vs_r_halo_disk_gas[s,0,:,selec] != 0)
    xplot = xlf[ind]
    yplot = r_vs_r_halo_disk_gas[s,0,ind,selec]
    errdn = r_vs_r_halo_disk_gas[s,1,ind,selec]
    errup = r_vs_r_halo_disk_gas[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='b', linestyle='dotted')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.2,interpolate=True)

    ax.plot(xlf,xlf+np.log10(0.017),color='DarkOrange',linewidth=1,label='Somerville+18')
    ax.plot(xlf,xlf+np.log10(0.019),color='DarkOrange',linewidth=1,linestyle='dotted')
    ax.plot(xlf,xlf+np.log10(0.016),color='DarkOrange',linewidth=1,linestyle='dotted')

    common.prepare_legend(ax, ['DarkOrange'], loc=2)

    common.savefig(outdir, fig, 'halo_galaxy_connection.pdf')
 

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'm_bh', 'rstar_disk', 'rstar_bulge', 'type', 
                           'specific_angular_momentum_disk_star', 'specific_angular_momentum_bulge_star',
                           'specific_angular_momentum_disk_gas', 'specific_angular_momentum_bulge_gas',
                           'specific_angular_momentum_disk_gas_atom', 'specific_angular_momentum_disk_gas_mol',
                           'lambda_subhalo', 'mvir_subhalo', 'mvir_hosthalo', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge','sfr_disk', 'sfr_burst', 'vvir_hosthalo', 
                           'rgas_disk', 'rgas_bulge')}

    # Loop over redshift and subvolumes

    sam_vs_sam_halo_disk     = np.zeros(shape = (len(zlist), 3, len(xlf), 2))
    sam_vs_sam_halo_gal      = np.zeros(shape = (len(zlist), 3, len(xlf), 2))
    sam_vs_sam_halo_disk_gas = np.zeros(shape = (len(zlist), 3, len(xlf), 2))
    sam_vs_sam_halo_bar      = np.zeros(shape = (len(zlist), 3, len(xlf), 2))

    m_vs_m_halo_disk     = np.zeros(shape = (len(zlist), 3, len(xmf), 2))
    m_vs_m_halo_gal      = np.zeros(shape = (len(zlist), 3, len(xmf), 2))
    m_vs_m_halo_disk_gas = np.zeros(shape = (len(zlist), 3, len(xmf), 2))
    m_vs_m_halo_bar      = np.zeros(shape = (len(zlist), 3, len(xmf), 2))

    r_vs_r_halo_disk     = np.zeros(shape = (len(zlist), 3, len(xlf), 2))
    r_vs_r_halo_gal      = np.zeros(shape = (len(zlist), 3, len(xlf), 2))
    r_vs_r_halo_disk_gas = np.zeros(shape = (len(zlist), 3, len(xlf), 2))
    r_vs_r_halo_bar      = np.zeros(shape = (len(zlist), 3, len(xlf), 2))

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        (lh, lj, lm, bt, ms, ssfr)  = prepare_data(hdf5_data, index, sam_vs_sam_halo_disk, sam_vs_sam_halo_gal,
                     sam_vs_sam_halo_disk_gas, sam_vs_sam_halo_bar, m_vs_m_halo_disk, m_vs_m_halo_gal, 
                     m_vs_m_halo_disk_gas, m_vs_m_halo_bar, r_vs_r_halo_disk, r_vs_r_halo_gal, 
                     r_vs_r_halo_disk_gas, r_vs_r_halo_bar)

    plot_specific_am_ratio(plt, outdir, obsdir, sam_vs_sam_halo_disk, sam_vs_sam_halo_gal, sam_vs_sam_halo_disk_gas, 
            sam_vs_sam_halo_bar, m_vs_m_halo_disk, m_vs_m_halo_gal,
            m_vs_m_halo_disk_gas, m_vs_m_halo_bar, r_vs_r_halo_disk, 
            r_vs_r_halo_gal, r_vs_r_halo_disk_gas, r_vs_r_halo_bar)

if __name__ == '__main__':
    main(*common.parse_args())
