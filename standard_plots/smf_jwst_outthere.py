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
mlow = 7
mupp = 12
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

def prepare_data(hdf5_data, seds_ap, index, hist_smf, hist_smf_kbright, hist_smf_sf, hist_smf_sf_kbright, hist_smf_q, hist_smf_q_kbright, z):

    (h0, volh, mdisk, mbulge, sfr_disk, sfr_bulge) = hdf5_data

    mag_ap_all = seds_ap[0] #apparent magnitudes with dust
    kband = 10
    kthresh = 24

    #From Furlong et al. (2015)
    if(z <= 2):
       ssfr_thresh = -2.0 + 0.5 * z
    else:
       ssfr_thresh = -1

    ind = np.where((mdisk+mbulge) > 0.0)
    mass = np.log10(mdisk[ind] + mbulge[ind]) - np.log10(float(h0))
    ssfr = (sfr_disk[ind] + sfr_bulge[ind]) / (mdisk[ind] + mbulge[ind]) #in Gyr^-1
 

    ind = np.where(mass > 0.0)
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf[index,:] = hist_smf[index,:] + H

    ind = np.where((mass > 0.0) & (mag_ap_all[kband,:] <= kthresh) & (mag_ap_all[kband,:] > 0))
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_kbright[index,:] = hist_smf_kbright[index,:] + H

    #select star-forming galaxies
    ind = np.where((mass > 0.0) & (ssfr > 10.0**ssfr_thresh))
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_sf[index,:] = hist_smf_sf[index,:] + H

    ind = np.where((mass > 0.0) & (mag_ap_all[kband,:] <= kthresh) & (mag_ap_all[kband,:] > 0) & (ssfr > 10.0**ssfr_thresh)) 
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_sf_kbright[index,:] = hist_smf_sf_kbright[index,:] + H

    #select quiescent galaxies
    ind = np.where((mass > 0.0) & (ssfr < 10.0**ssfr_thresh))
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_q[index,:] = hist_smf_q[index,:] + H

    ind = np.where((mass > 0.0) & (mag_ap_all[kband,:] <= kthresh) & (mag_ap_all[kband,:] > 0) & (ssfr < 10.0**ssfr_thresh)) 
    H, _ = np.histogram(mass[ind],bins=np.append(mbins,mupp))
    hist_smf_q_kbright[index,:] = hist_smf_q_kbright[index,:] + H

    if volh > 0:
        vol = volh/pow(h0,3.)  # In Mpc^3
        hist_smf[index,:]  = hist_smf[index,:]/vol/dm
        hist_smf_kbright[index,:]  = hist_smf_kbright[index,:]/vol/dm
        hist_smf_sf[index,:]  = hist_smf_sf[index,:]/vol/dm
        hist_smf_sf_kbright[index,:]  = hist_smf_sf_kbright[index,:]/vol/dm
        hist_smf_q[index,:]  = hist_smf_q[index,:]/vol/dm
        hist_smf_q_kbright[index,:]  = hist_smf_q_kbright[index,:]/vol/dm

    return mass


def plot_smf_jwst(obsdir, outdir, plt, zlist, hist_smf, hist_smf_kbright, hist_smf_sf, hist_smf_sf_kbright, hist_smf_q, hist_smf_q_kbright):

    fig = plt.figure(figsize=(11.7,11.7))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm M_{\\star}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 8, 13, -6, -1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (331, 332, 333, 334, 335, 336, 337, 338, 339)
    indeces = (0, 1, 2, 3, 4, 5, 6, 7, 8)

    for subplot, idx, z in zip(subplots, indeces, zlist):

        ax = fig.add_subplot(subplot)
        if(idx == 0 or idx == 3 or idx == 6):
           ytitle = ytit
        else:
           ytitle = ' '
        
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(z))

        y = hist_smf[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'k', linestyle='solid', label='all' if idx == 0 else None)
        y = hist_smf_kbright[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'k', linestyle='dotted', label='K<24' if idx == 0 else None)
        y = hist_smf_sf[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'b', linestyle='solid', label='SF' if idx == 0 else None)
        y = hist_smf_sf_kbright[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'b', linestyle='dotted', label='SF K<24' if idx == 0 else None)
        y = hist_smf_q[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'r', linestyle='solid', label='Quiesc' if idx == 0 else None)
        y = hist_smf_q_kbright[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'r', linestyle='dotted', label='Quiesc K<24' if idx == 0 else None)

        if(idx == 0):
           common.prepare_legend(ax, ['k','k','b','b','r','r'], loc = 'upper right')

    common.savefig(outdir, fig, 'stellarmf_z_jwst_outthere.pdf')


def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist = (0.9, 1.59696, 2.00392, 2.47464723643932, 3.01916, 3.50099697082904, 3.95972, 4.465197621546, 5.02220991014863)

    plt = common.load_matplotlib()

    # Histograms
    hist_smf         = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_kbright = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_sf         = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_sf_kbright = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_q          = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_q_kbright  = np.zeros(shape = (len(zlist), len(mbins)))

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'sfr_disk', 'sfr_burst')}

    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5"
    fields_sed_ap = {'SED/ap_dust': ('total','disk'),}

    #Bands information:
    #(0): "FUV_GALEX", "NUV_GALEX", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS",
    #(6): "z_SDSS", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1_WISE",
    #(12): "I1_Spitzer", "I2_Spitzer", "W2_WISE", "I3_Spitzer", "I4_Spitzer",
    #(17): "W3_WISE", "W4_WISE", "P70_Herschel", "P100_Herschel",
    #(21): "P160_Herschel", "S250_Herschel", "S350_Herschel", "S450_JCMT",
    #(25): "S500_Herschel", "S850_JCMT", "Band9_ALMA", "Band8_ALMA",
    #(29): "Band7_ALMA", "Band6_ALMA", "Band5_ALMA", "Band4_ALMA"

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        seds_ap = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed_ap, subvols, file_hdf5_sed)

        mass = prepare_data(hdf5_data, seds_ap, index, hist_smf, hist_smf_kbright, hist_smf_sf, hist_smf_sf_kbright, hist_smf_q, hist_smf_q_kbright, zlist[index])
        h0 = hdf5_data[0]

    # Take logs
    def take_log(array):
        ind = np.where(array > 0.)
        array[ind] = np.log10(array[ind])

    take_log(hist_smf)
    take_log(hist_smf_kbright)
    take_log(hist_smf_sf)
    take_log(hist_smf_sf_kbright)
    take_log(hist_smf_q)
    take_log(hist_smf_q_kbright)

    for index, snapshot in enumerate(redshift_table[zlist]):
        print("#SMF at z=%s and snapshot=%s" % (str(zlist[index]), str(redshift_table[zlist[index]])))
        print("#phi in units of log10(Mpc^-3 dex^-1)")
        print("#log10(Mstar/Msun) phi_all phi_all_kbright phi_sf phi_sf_kbright phi_quiesc phi_quiesc_kbright")
        for a,b,c,d,e,f,g in zip(xmf, hist_smf[index,:], hist_smf_kbright[index,:], hist_smf_sf[index,:], hist_smf_sf_kbright[index,:], hist_smf_q[index,:], hist_smf_q_kbright[index,:]):
            print(a,b,c,d,e,f,g)

    plot_smf_jwst(obsdir, outdir, plt, zlist, hist_smf, hist_smf_kbright, hist_smf_sf, hist_smf_sf_kbright, hist_smf_q, hist_smf_q_kbright)
 
if __name__ == '__main__':
    main(*common.parse_args())
