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
PI       = 3.1416

mlow = 6.5
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

r_ap = 0.1 #in units of Mpc

rlow = 0
rupp = 1.0
dr = 0.005
rbins = np.arange(rlow,rupp,dr)
xrf = rbins + dr/2.0

use_r50_aperture = True
r50_aperture = 10.0

def prepare_data(hdf5_data, index, rcomb, rcomb_icl, bs_error):

    (h0, _, mdisk, mbulge, rdisk, rbulge, typeg, mstellar_halo, cnfw, vvir, mvir) = hdf5_data

    def correct_h(x,h0):
        return x/h0

    mstars_tot = correct_h((mdisk+mbulge), h0)
    rdisk = correct_h(rdisk,h0)
    rbulge = correct_h(rbulge,h0)
    mdisk = correct_h(mdisk,h0)
    mbulge = correct_h(mbulge,h0)
    mstellar_halo = correct_h(mstellar_halo,h0)

    bin_it   = functools.partial(us.wmedians, xbins=xmf)

    #compute a stellar mass weighted size for the galaxy
    rcombined = (mdisk * rdisk  + mbulge * rbulge) / (mdisk + mbulge)
    
    #calculate size-mass relation for rcombined
    ind = np.where(mdisk+mbulge > 0)
    rcomb[index,:] = bin_it(x=np.log10(mstars_tot[ind]),
                            y=np.log10(rcombined[ind]*MpcToKpc))

    #model inclusion of ICL in sizes
    rnew = np.zeros(shape = len(rcombined))
    def enclosed_mass_icl(x, rho_icl, c, rvir):
        return 4.0 * PI * rho_icl * (c * rvir)**3.0 * (1.0 / (1.0 + c * x) - 1 + np.log(1.0 + c * x))

    rvir  = G * mvir / pow(vvir,2.0) / h0
    concen_factor = np.log(1.0 + cnfw) - cnfw / (1.0 + cnfw)

    #find centrals with stellar halos
    ind = np.where((typeg == 0) & (mstellar_halo > 0))
    cnfw_in = cnfw[ind]
    rvir_in = rvir[ind]
    rho_icl = mstellar_halo[ind]/ (4.0 * PI * (cnfw_in * rvir_in)**3.0) / concen_factor[ind]
    rcombined_in = rcombined[ind] #in Mpc

    #compute aperture at which ther half-mass of the ICL will be computed.
    x = r_ap / rvir_in #fixed aperture

    if(use_r50_aperture == True): #compute aperture based on r50
       x = r50_aperture * rcombined_in / rvir_in

    micl_apperture = np.zeros(shape = len(rho_icl))
    for i in range(0, len(rho_icl)):
        micl_apperture[i] = enclosed_mass_icl(x[i], rho_icl[i], cnfw_in[i], rvir_in[i])
    micl_apperture_50 = micl_apperture * 0.5
    
    #use a linear interpolation to find r50, between log10(mass) and r/rvir.
    def find_r50_icl(xrf, yrf, mass):
        find_nearest = np.argsort(abs(yrf-mass))
        m = (xrf[find_nearest[1]] - xrf[find_nearest[0]]) / (yrf[find_nearest[1]] - yrf[find_nearest[0]])
        a = xrf[find_nearest[1]] - m * yrf[find_nearest[1]]
        return m * mass + a

    #compute r50_icl for each galaxy
    r50_icl = np.zeros(shape = len(rho_icl))
    
    for i in range(0, len(rho_icl)):
        r50_icl[i] = find_r50_icl(xrf, np.log10(enclosed_mass_icl(xrf, rho_icl[i], cnfw_in[i], rvir_in[i])), np.log10(micl_apperture_50[i]))
    
    
    #compute a new size doing a stellar mass weighted size
    rnew[ind] = (rcombined[ind] * mstars_tot[ind] + r50_icl * micl_apperture) / (mstars_tot[ind] + micl_apperture)
    
    #for galaxies that have not been assigned an rnew, we use rcombined (those include satellites and centrals without stellar halos)
    ind = np.where(rnew == 0)
    rnew[ind] = rcombined[ind]
   
    #compute the size-mass relation now including the icl.
    ind = np.where(mdisk+mbulge > 0)
    rcomb_icl[index,:] = bin_it(x=np.log10(mstars_tot[ind]),
                            y=np.log10(rnew[ind]*MpcToKpc))
    # compute bootstrapped error on updated median sizes
    bs_error[index] = us.bootstrap_error(x=np.log10(mstars_tot[ind]),
                            y=np.log10(rnew[ind]*MpcToKpc), xbins=xmf)
   

def plot_sizes_combined(plt, outdir, rcomb, rcomb_icl, bs_error, obsdir):
    fig = plt.figure(figsize=(5,4.5))
    
    # Total ##################################
    xtit="$\\rm log_{10} (\\rm M_{\\rm stars, total}/M_{\odot})$"
    ytit="$\\rm log_{10} (\\rm r_{\\rm eff, comb}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.5, 2

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    #Predicted size-mass for disks
    ind = np.where(rcomb[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = rcomb[0,0,ind] # z=0 and median only
    errdn = rcomb[0,1,ind]
    errup = rcomb[0,2,ind]
    #ax.errorbar(xplot,yplot[0],yerr=bs_error[:,0,ind][0], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o',label="Shark galaxies")
    ax.plot(xplot, yplot[0],  ls='None', mfc='None', mec='k',marker='o',label="Shark galaxies")

    ind = np.where(rcomb_icl[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = rcomb_icl[0,0,ind] # z=0 and median only
    #errdn = rcomb_icl[0,1,ind]
    #errup = rcomb_icl[0,2,ind]
    ax.errorbar(xplot,yplot[0],yerr=bs_error[0,ind], ls='None', mfc='None', ecolor = 'r', mec='r',marker='s',label="Shark galaxies+icl")

    # load semi-major sizes
    obs = common.load_observation(obsdir,'SizeMass/GAMA_H-band_dlogM_0.25_reff.txt', [0,1,2,3])
    xobs = obs[0]
    yobs = obs[1]
    errobs = obs[3]
    ax.plot(xobs, yobs, color = 'teal', label = 'GAMA H-band')
    ax.fill_between(xobs, yobs-errobs, yobs+errobs, color = 'teal', alpha=0.4)
    
    common.prepare_legend(ax, ['k','k','k'], loc=2)
    common.savefig(outdir, fig, 'sizes_combined_ICL.pdf')

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'rstar_disk', 'rstar_bulge', 'type', 'mstellar_halo', 'cnfw_subhalo', 'vvir_hosthalo', 'mvir_hosthalo')}

    # Loop over redshift and subvolumes
    rcomb = np.zeros(shape = (len(zlist), 3, len(xmf)))
    rcomb_icl = np.zeros(shape = (len(zlist), 3, len(xmf)))
    bs_error = np.zeros(shape = (len(zlist), len(xmf)))

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        prepare_data(hdf5_data, index, rcomb, rcomb_icl, bs_error)
 
    plot_sizes_combined(plt, outdir, rcomb, rcomb_icl, bs_error, obsdir)
if __name__ == '__main__':
    main(*common.parse_args())
