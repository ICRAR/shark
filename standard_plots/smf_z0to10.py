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
"""SMF plots"""

import collections
import functools
import logging
import math

import numpy as np

import common
import utilities_statistics as us


observation = collections.namedtuple('observation', 'label x y yerrup yerrdn err_absolute')

logger = logging.getLogger(__name__)

##################################
# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
MpcToKpc = 1e3

##################################
# Mass function initialization
mlow = 5
mupp = 14
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0
imf   = 'cha'

mlow2 = 5
mupp2 = 14
dm2 = 0.3
mbins2 = np.arange(mlow2,mupp2,dm2)
xmf2 = mbins2 + dm2/2.0

ssfrlow = -6
ssfrupp = 4
dssfr = 0.2
ssfrbins = np.arange(ssfrlow,ssfrupp,dssfr)
xssfr    = ssfrbins + dssfr/2.0

sfrlow = -3
sfrupp = 1.5
dsfr = 0.2
sfrbins = np.arange(sfrlow,sfrupp,dsfr)
xsfr    = sfrbins + dsfr/2.0


def load_smf_observations(obsdir, h0):
    # Thorne et al. (2021)
    def add_thorne21_data(zappend, file_read='z0.4'):
        lm, pD, dn, du = np.loadtxt(obsdir+'/mf/SMF/Thorne21/SMFvals_'+file_read+'.csv', delimiter=',', skiprows=1, usecols = [0,1,2,3], unpack = True)
        hobd = 0.7
        pDlog = np.log10(pD[::3]) +  3.0 * np.log10(hobs/h0)
        dnlog = np.log10(pD[::3]) - np.log10(dn[::3])
        dulog = np.log10(du[::3]) - np.log10(pD[::3])
        lm = lm[::3] -  2.0 * np.log10(hobs/h0)
        zappend.append((observation("Thorne+2021", lm, pDlog, dnlog, dulog, err_absolute=False), 's'))

    # Weaver et al. (2022; COSMOS2020)
    def add_weaver22_data(zappend, file_read='0.2z0.5'):
        lm, pD, dn, du = np.loadtxt(obsdir+'/mf/SMF/COSMOS2020/SMF_Farmer_v2.1_' + file_read + '_total.txt', delimiter=' ', skiprows=0, usecols = [0,2,3,4], unpack = True)
        hobd = 0.7
        #find cases with negative number densities and change those with a very small value 
        ind = np.where(dn < 0)
        dn[ind] = 1e-10
        pDlog = np.log10(pD) +  3.0 * np.log10(hobs/h0)
        dnlog = np.log10(pD) - np.log10(dn)
        dulog = np.log10(du) - np.log10(pD)
        lm = lm -  2.0 * np.log10(hobs/h0)
        zappend.append((observation("Weaver+2023", lm, pDlog, dnlog, dulog, err_absolute=False), '*'))


    # Driver al. (2022, z=0). Chabrier IMF
    z0obs = []
    lm, p, dp = common.load_observation(obsdir, 'mf/SMF/GAMAIV_Driver22.dat', [0,1,2])
    hobs = 0.7
    xobs = lm + 2.0 * np.log10(hobs/h0)
    yobs = p - 3.0 * np.log10(hobs/h0)
    z0obs.append((observation("Driver+2022", xobs, yobs, dp, dp, err_absolute=False), 'o'))

    lm, p, dpdn, dpup = common.load_observation(obsdir, 'mf/SMF/SMF_Bernardi2013_SerExp.data', [0,1,2,3])
    xobs = lm + 2.0 * np.log10(hobs/h0)
    yobs = np.log10(p) - 3.0 * np.log10(hobs/h0)
    ydn = np.log10(p) - np.log10(p-dpdn)
    yup = np.log10(p+dpup) - np.log10(p) 
    z0obs.append((observation("Bernardi+2013", xobs, yobs, ydn, yup, err_absolute=False), 's'))


    lm, p, dpdn, dpup = common.load_observation(obsdir, 'mf/SMF/SMF_Li2009.dat', [0,1,2,3])
    xobs = lm - 2.0 * np.log10(hobs) + 2.0 * np.log10(hobs/h0)
    yobs = p + 3.0 * np.log10(hobs) - 3.0 * np.log10(hobs/h0)
    z0obs.append((observation("Li&White+2009", xobs, yobs, abs(dpdn), dpup, err_absolute=False), 'd'))


    # Moustakas (Chabrier IMF), ['Moustakas+2013, several redshifts']
    zdnM13, lmM13, pM13, dp_dn_M13, dp_up_M13 = common.load_observation(obsdir, 'mf/SMF/SMF_Moustakas2013.dat', [0,3,5,6,7])
    xobsM13 = lmM13 + 2.0 * np.log10(hobs/h0)

    yobsM13 = np.full(xobsM13.shape, -999.) - 3.0 * np.log10(hobs/h0)
    lerrM13 = np.full(xobsM13.shape, -999.)
    herrM13 = np.full(xobsM13.shape, 999.)
    indx = np.where( pM13 < 1)
    yobsM13[indx] = (pM13[indx])
    indx = np.where( dp_dn_M13 > 0)
    lerrM13[indx]  = dp_dn_M13[indx] 
    indx = np.where( dp_up_M13 > 0)
    herrM13[indx]  = dp_up_M13[indx]

    # Muzzin (Kroupa IMF), ['Moustakas+2013, several redshifts']
    zdnMu13,zupMu13,lmMu13,pMu13,dp_dn_Mu13,dp_up_Mu13 = common.load_observation(obsdir, 'mf/SMF/SMF_Muzzin2013.dat', [0,1,2,4,5,5])
    # -0.09 corresponds to the IMF correction
    xobsMu13 = lmMu13 - 0.09 + 2.0 * np.log10(hobs/h0) 
    yobsMu13 = np.full(xobsMu13.shape, -999.) - 3.0 * np.log10(hobs/h0)
    lerrMu13 = np.full(xobsMu13.shape, -999.)
    herrMu13 = np.full(xobsMu13.shape, 999.)
    indx = np.where( pMu13 < 1)
    yobsMu13[indx] = (pMu13[indx])
    indx = np.where( dp_dn_Mu13 > 0)
    lerrMu13[indx]  = dp_dn_Mu13[indx] 
    indx = np.where( dp_up_Mu13 > 0)
    herrMu13[indx]  = dp_up_Mu13[indx]

    # z0.5 obs
    z05obs = []
    #in_redshift = np.where(zdnM13 == 0.4)
    #z05obs.append((observation("Moustakas+2013", xobsM13[in_redshift], yobsM13[in_redshift], lerrM13[in_redshift], herrM13[in_redshift], err_absolute=False), 'o'))
    in_redshift = np.where(zdnMu13 == 0.5)
    z05obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z05obs, file_read='z0.51')
    add_weaver22_data(z05obs, file_read='0.2z0.5')


    # z1 obs
    z1obs = []
    #in_redshift = np.where(zdnM13 == 0.8)
    #z1obs.append((observation("Moustakas+2013", xobsM13[in_redshift], yobsM13[in_redshift], lerrM13[in_redshift], herrM13[in_redshift], err_absolute=False), 'o'))
    in_redshift = np.where(zdnMu13 == 1)
    z1obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z1obs, file_read='z1.1')
    add_weaver22_data(z1obs, file_read='0.8z1.1')

    #z2 obs
    z2obs = []
    in_redshift = np.where(zupMu13 == 2.5)
    z2obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    #in_redshift = np.where(zdnS12 == 1.8)
    #z2obs.append((observation("Santini+2012", xobsS12[in_redshift], yobsS12[in_redshift], lerrS12[in_redshift], herrS12[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z2obs, file_read='z2')
    add_weaver22_data(z2obs, file_read='2.0z2.5')

    # z3 obs
    z3obs = []
    in_redshift = np.where(zupMu13 == 3.0)
    z3obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    #in_redshift = np.where(zdnS12 == 2.5)
    #z3obs.append((observation("Santini+2012", xobsS12[in_redshift], yobsS12[in_redshift], lerrS12[in_redshift], herrS12[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z3obs, file_read='z3')
    add_weaver22_data(z3obs, file_read='3.0z3.5')

    # z4 obs
    z4obs = []
    in_redshift = np.where(zupMu13 == 4.0)
    z4obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    #in_redshift = np.where(zdnS12 == 3.5)
    #z4obs.append((observation("Santini+2012", xobsS12[in_redshift], yobsS12[in_redshift], lerrS12[in_redshift], herrS12[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z4obs, file_read='z4')
    add_weaver22_data(z4obs, file_read='3.5z4.5')


    # z6 obs
    z6obs = []
    add_weaver22_data(z6obs, file_read='5.5z6.5')

    # z8 obs
    z8obs = []
    add_weaver22_data(z7obs, file_read='6.5z7.5')

    # z10 obs
    z10obs = []

    return (z0obs, z05obs, z1obs, z2obs, z3obs, z4obs, z6obs, z8obs, z10obs)

def plot_stellarmf_z(plt, outdir, obsdir, h0, plotz, hist_smf, hist_smf_cen, hist_smf_sat, hist_smf_offset, hist_smf_30kpc):

    (z0obs, z05obs, z1obs, z2obs, z3obs, z4obs, z6obs, z8obs, z10obs) = load_smf_observations(obsdir, h0)
    
    PlotLagos18 = True
    def plot_lagos18_smf(ax, z):
        sm, z0, z0p5, z1, z2, z3, z4 = common.load_observation(obsdir, 'Models/SharkVariations/SMF_Lagos18.dat', [0,1,2,3,4,5,6])
        if z == 0:
           y = z0
        elif z == 1:
           y = z0p5
        elif z == 2:
           y = z1
        elif z == 3:
           y = z2
        elif z == 4:
           y = z3
        elif z == 5:
           y = z4
        ax.plot(sm, y, linestyle='dashed', color='black',label='Shark v1.1 (L18)' if z == 0 else None)

    
    fig = plt.figure(figsize=(9.7,11.7))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm M_{\\star}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 8, 13, -6, -1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (331, 332, 333, 334, 335, 336, 337, 338, 339)
    indeces = (0, 1, 2, 3, 4, 5, 6, 7, 8)
    zs = (0, 0.5, 1, 2, 3, 4, 6, 8, 10)
    observations = (z0obs, z05obs, z1obs, z2obs, z3obs, z4obs, z6obs, z8obs, z10obs)

    for subplot, idx, z, obs_and_markers in zip(subplots, indeces, zs, observations):

        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(z))

        # Observations
        for obs, marker in obs_and_markers:
            common.errorbars(ax, obs.x, obs.y, obs.yerrdn, obs.yerrup, 'grey',
                             marker, err_absolute=obs.err_absolute, label=obs.label, markersize=4)

        # Predicted SMF
        if plotz[idx]:
            y = hist_smf[idx,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'r', label='Shark v2.0' if idx == 0 else None)
            #print('#stellar mass function redshift', z)
            #for a,b in zip(xmf[ind],y[ind]):
            #    print(a,b)
            #y = hist_smf_cen[idx,:]
            #ind = np.where(y < 0.)
            #ax.plot(xmf[ind],y[ind],'b', linestyle='dotted', label ='centrals' if idx == 0 else None)
            #y = hist_smf_sat[idx,:]
            #ind = np.where(y < 0.)
            #ax.plot(xmf[ind],y[ind],'g', linestyle='dotted', label ='Shark v2.0 (sats)' if idx == 0 else None)
            if z < 1:
                y = hist_smf_30kpc[idx,:]
                ind = np.where(y < 0.)
                ax.plot(xmf[ind],y[ind],'r', linestyle='dotted', linewidth=1, label ='Shark v2.0 (30kpc)'  if idx == 0 else None)
            if z >= 1:
                y = hist_smf_offset[idx,:]
                ind = np.where(y < 0.)
                ax.plot(xmf[ind],y[ind],'r', linestyle='dashdot', linewidth=2, label ='0.3dex error')
            if PlotLagos18 == True:
               if(z <= 4):
                  plot_lagos18_smf(ax, idx)
       

        colors = []
        if idx == 0:
            colors = ['r','r','k']
        #elif idx == 1:
        #    colors += ['r']
        elif idx > 1:
            colors = ['r']
        colors += ['grey', 'grey','grey']

        common.prepare_legend(ax, colors)

    common.savefig(outdir, fig, 'stellarmf_z0to10.pdf')


def prepare_data(hdf5_data, index, hist_smf, hist_smf_offset, hist_smf_cen, hist_smf_sat, 
                 hist_smf_30kpc, plotz, hist_smf_err, hist_smf_comp):

    (h0, volh, sfr_disk, sfr_burst, mdisk, mbulge, rstar_disk, mBH, mHI, mH2, 
     mgas_disk, mHI_bulge, mH2_bulge, mgas_bulge, mgas_metals_disk, mgas_metals_bulge, 
     mstars_metals_disk, mstars_metals_bulge, typeg, mvir_hosthalo, rstar_bulge, 
     mbulge_mergers, mbulge_diskins, mbulge_mergers_assembly, mbulge_diskins_assembly) = hdf5_data

    zstar = (mstars_metals_disk + mstars_metals_bulge) / (mdisk + mbulge)
    rcomb = (rstar_disk * mdisk + rstar_bulge * mbulge) / (mdisk + mbulge) / h0 * 1e3

    mass          = np.zeros(shape = len(mdisk))
    mass_30kpc    = np.zeros(shape = len(mdisk))
    massd_30kpc   = np.zeros(shape = len(mdisk))
    massb_30kpc   = np.zeros(shape = len(mdisk))


    ind = np.where((np.isnan(mdisk+mbulge) == True) | (np.isinf(mdisk+mbulge) == True))
    if (len(mdisk[ind]) > 0):
         print("Number of galaxies with a stellar mass of NaN:", len(mdisk[ind]))
    else:
         print("All galaxies have well defined stellar mass")

    ind = np.where((np.isnan(mBH) == True) | (np.isinf(mBH) == True))
    if (len(mdisk[ind]) > 0):
         print("Number of galaxies with a BH mass of NaN:", len(mdisk[ind]))
    else:
         print("All galaxies have well defined BH mass")


    #select massive centrals

    ind = np.where((mdisk+mbulge) > 0.0)
    mass[ind] = np.log10(mdisk[ind] + mbulge[ind]) - np.log10(float(h0))
    logger.debug('number of galaxies with mstars>0 and max mass: %d, %d', len(mass[ind]), max(mass[ind]))
    
    H, _ = np.histogram(mass,bins=np.append(mbins,mupp))
    hist_smf[index,:] = hist_smf[index,:] + H
    ran_err = np.random.normal(0.0, 0.3, len(mass))
    mass_err = mass + ran_err
    H, _ = np.histogram(mass_err,bins=np.append(mbins,mupp))
    hist_smf_offset[index,:] = hist_smf_offset[index,:] + H

    #Calculate the stellar mass contained in 30pkpc, assuming an exponential profile for the disk and a Plummer profile for the bulge.
    ind = np.where((mdisk > 0.0)  & (rstar_disk > 0))
    massd_30kpc[ind] = mdisk[ind] * (1.0 - (1.0 + 30.0/(rstar_disk[ind]/1.67/h0 * MpcToKpc)) * np.exp(-30.0/(rstar_disk[ind]/1.67/h0 * MpcToKpc)))
    ind = np.where((mbulge > 0.0)  & (rstar_bulge > 0))
    massb_30kpc[ind] = mbulge[ind] * pow(30.0, 3.0) / pow((pow(30.0, 2.0) + pow(rstar_bulge[ind]/1.3/h0 * MpcToKpc, 2.0)), 3.0/2.0)

    ind = np.where((massd_30kpc + massb_30kpc) > 0)
    mass_30kpc[ind] = np.log10(massd_30kpc[ind] + massb_30kpc[ind]) - np.log10(float(h0))
    H, _ = np.histogram(mass_30kpc,bins=np.append(mbins,mupp))
    hist_smf_30kpc[index,:] = hist_smf_30kpc[index,:] + H

    #stellar mass functions separated into centrals and satellites
    ind = np.where(typeg == 0)
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_cen[index,:] = hist_smf_cen[index,:] + H
    ind = np.where(typeg > 0)
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_sat[index,:] = hist_smf_sat[index,:] + H

    #stellar mass functions by galaxy components (disks, bulges, bulges by mergers, bulges by disk ins)
    ind = np.where(mdisk > 0)
    H, _ = np.histogram(np.log10(mdisk[ind]),bins=np.append(mbins,mupp))
    hist_smf_comp[index,0,:] = hist_smf_comp[index,0,:] + H
    ind = np.where(mbulge > 0)
    H, _ = np.histogram(np.log10(mbulge[ind]),bins=np.append(mbins,mupp))
    hist_smf_comp[index,1,:] = hist_smf_comp[index,1,:] + H
    ind = np.where(mbulge_mergers > 0)
    H, _ = np.histogram(np.log10(mbulge_mergers[ind]),bins=np.append(mbins,mupp))
    hist_smf_comp[index,2,:] = hist_smf_comp[index,2,:] + H
    ind = np.where(mbulge_diskins > 0)
    H, _ = np.histogram(np.log10(mbulge_diskins[ind]),bins=np.append(mbins,mupp))
    hist_smf_comp[index,3,:] = hist_smf_comp[index,3,:] + H

    if volh > 0:
        vol = volh/pow(h0,3.)  # In Mpc^3
        hist_smf_err[index,:]  = (hist_smf[index,:] - np.sqrt(hist_smf[index,:]))/vol/dm

        hist_smf[index,:]  = hist_smf[index,:]/vol/dm
        hist_smf_comp[index,:]  = hist_smf_comp[index,:]/vol/dm
        hist_smf_30kpc[index,:] = hist_smf_30kpc[index,:]/vol/dm
        hist_smf_offset[index,:] = hist_smf_offset[index,:]/vol/dm
        hist_smf_cen[index,:]  = hist_smf_cen[index,:]/vol/dm
        hist_smf_sat[index,:]  = hist_smf_sat[index,:]/vol/dm

        plotz[index]     = True
    else:
        plotz[index]     = False

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist = (0, 0.5, 1, 2, 3, 4, 6, 8, 10)
    #zlist = (0.005, 0.2, 0.5 , 0.8 , 1.1 , 1.5 , 2.2 , 2.9 , 3.9, 5.1)

    plt = common.load_matplotlib()

    # Histograms
    hist_smf       = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_30kpc = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_offset   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_cen   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_sat   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_err   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_comp  = np.zeros(shape = (len(zlist), 4, len(mbins)))

    plotz = np.empty(shape=(len(zlist)), dtype=np.bool_)

    fields = {'galaxies': ('sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge',
                           'mgas_metals_disk', 'mgas_metals_bulge',
                           'mstars_metals_disk', 'mstars_metals_bulge', 'type', 
                           'mvir_hosthalo', 'rstar_bulge', 'mstars_burst_mergers', 
                           'mstars_burst_diskinstabilities', 'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        prepare_data(hdf5_data, index, hist_smf, hist_smf_offset, hist_smf_cen,
                     hist_smf_sat, hist_smf_30kpc, plotz, hist_smf_err, hist_smf_comp)

        h0 = hdf5_data[0]
        volh = hdf5_data[1]

    # This should be the same in all HDF5 files

    # Take logs
    def take_log(array):
        ind = np.where(array > 0.)
        array[ind] = np.log10(array[ind])

    take_log(hist_smf)
    take_log(hist_smf_comp)
    take_log(hist_smf_30kpc)
    take_log(hist_smf_cen)
    take_log(hist_smf_sat)
    take_log(hist_smf_offset)

    plot_stellarmf_z(plt, outdir, obsdir, h0, plotz, hist_smf, hist_smf_cen, hist_smf_sat, hist_smf_offset, hist_smf_30kpc)

if __name__ == '__main__':
    main(*common.parse_args())
