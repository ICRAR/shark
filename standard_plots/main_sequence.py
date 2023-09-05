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
dm = 0.125
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0
imf   = 'cha'

mlow2 = 5
mupp2 = 14
dm2 = 0.3
mbins2 = np.arange(mlow2,mupp2,dm2)
xmf2 = mbins2 + dm2/2.0

mlow3 = 5
mupp3 = 14
dm3 = 0.2
mbins3 = np.arange(mlow3,mupp3,dm3)
xmf3 = mbins3 + dm3/2.0

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

def plot_sfr_mstars_z0(plt, outdir, obsdir, h0, sfr_seq, mainseqsf):

    bin_it = functools.partial(us.wmedians, xbins=xmf, nmin=50)

    fig = plt.figure(figsize=(5,5))
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(\\rm SFR/M_{\odot} yr^{-1})$"

    xmin, xmax, ymin, ymax = 8, 12.5, -3.2, 3.3
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    ind = np.where(sfr_seq[1,:] < -3)
    sfr_seq[1,ind] = -3
    #predicted relation
    ind = np.where((sfr_seq[0,:] > 7) & (sfr_seq[0,:] < 13) & (sfr_seq[1,:] >-10) & (sfr_seq[1,:] < 10))
    xdata = sfr_seq[0,ind]
    ydata = sfr_seq[1,ind]
    us.density_contour(ax, xdata[0], ydata[0], 30, 30, cmap = 'cividis') #, **contour_kwargs)

    ind = np.where(sfr_seq[0,:] > 0)
    toplot = bin_it(x=sfr_seq[0,ind], y=sfr_seq[1,ind])
    ind = np.where(toplot[0,:] != 0)
    yp = toplot[0,ind] 
    ax.plot(xmf[ind], yp[0],color='k',linestyle='solid', linewidth = 2, label="Shark v2.0")

    xm, ym = common.load_observation(obsdir, 'Models/SharkVariations/SFRMstars_Lagos18.dat', [0,1])
    ax.plot(xm, ym, linestyle='dashed', linewidth = 2, color='black',label='Shark v1.1 (L18)')

    #print("will print median MS")
    #for a,b in zip(xmf, toplot[0]):
    #    print(a,b)

    #SFR relation z=0
    lm, SFR = common.load_observation(obsdir, 'SFR/Brinchmann04.dat', (0, 1))
    hobs = 0.7
    #add cosmology correction plus IMF correction that goes into the stellar mass.
    corr_cos = np.log10(pow(hobs,2)/pow(h0,2)) - 0.09
    #apply correction to both stellar mass and SFRs.
    ax.plot(lm[0:35] + corr_cos, SFR[0:35] + corr_cos, color='SandyBrown', linewidth = 4, linestyle='dashed', label='Brinchmann+04')
    #ax.plot(lm[36:70] + corr_cos, SFR[36:70] + corr_cos, color='PaleVioletRed',linewidth = 5, linestyle='dotted')
    #ax.plot(lm[71:len(SFR)] + corr_cos, SFR[71:len(SFR)] + corr_cos, color='PaleVioletRed',linewidth = 5, linestyle='dotted')

    xdataD16 = [9.3, 10.6]
    ydataD16 = [-0.39, 0.477]
    ax.plot(xdataD16,ydataD16, color='Crimson',linestyle='dashdot',linewidth = 4, label='Davies+16')

    #GAMA data at z<0.06
    #CATAID StellarMass_bestfit StellarMass_50 StellarMass_16 StellarMass_84 SFR_bestfit SFR_50 SFR_16 SFR_84 Zgas_bestfit Zgas_50 Zgas_16 Zgas_84 DustMass_bestfit DustMass_50 DustMass_16 DustMass_84 DustLum_50 DustLum_16 DustLum_84 uberID redshift
    ms_gama, sfr_gama = common.load_observation(obsdir, 'GAMA/ProSpect_Claudia.txt', [2,6])
    ind = np.where(sfr_gama < 1e-3)
    sfr_gama[ind] = 1e-3
    #ax.hexbin(np.log10(ms_gama), np.log10(sfr_gama), gridsize=(20,20), mincnt=5) #, cmap = 'plasma') #, **contour_kwargs)
    us.density_contour_reduced(ax, np.log10(ms_gama), np.log10(sfr_gama), 25, 25) #, **contour_kwargs)

    toplot = bin_it(x=np.log10(ms_gama), y=np.log10(sfr_gama))
    ind = np.where(toplot[0,:] != 0)
    yp = toplot[0,ind]
    yup = toplot[2,ind]
    ydn = toplot[1,ind]
    ax.plot(xmf[ind], yp[0],color='Maroon',linestyle='dashed', linewidth = 5, label="Bellstedt+20")
#ax.plot(xmf[ind], yp[0]+yup[0],color='PaleVioletRed',linestyle='dotted', linewidth = 5)
    #ax.plot(xmf[ind], yp[0]-ydn[0],color='PaleVioletRed',linestyle='dotted', linewidth = 5)

    # individual massive galaxies from Terrazas+17
    ms, sfr, upperlimflag = common.load_observation(obsdir, 'BHs/MBH_host_gals_Terrazas17.dat', [0,1,2])
    ind = np.where(ms > 11.3)
    ax.errorbar(ms[ind], sfr[ind], xerr=0.2, yerr=0.3, ls='None', mfc='None', ecolor = 'r', mec='r',marker='s',label="Terrazas+17")
    ind = np.where((upperlimflag == 1) & (ms > 11.3))
    for a,b in zip (ms[ind], sfr[ind]):
        ax.arrow(a, b, 0, -0.3, head_width=0.05, head_length=0.1, fc='r', ec='r')

    # Legend
    common.prepare_legend(ax, ['k','k','SandyBrown','Crimson','Maroon','r'], loc=2)
    plt.tight_layout()
    common.savefig(outdir, fig, 'SFR_Mstars_z0.pdf')


def prepare_data(hdf5_data, index, mainseqsf):

    (h0, volh, sfr_disk, sfr_burst, mdisk, mbulge, rstar_disk, mBH, mHI, mH2, 
     mgas_disk, mHI_bulge, mH2_bulge, mgas_bulge, mgas_metals_disk, mgas_metals_bulge, 
     mstars_metals_disk, mstars_metals_bulge, typeg, mvir_hosthalo, rstar_bulge, 
     mbulge_mergers, mbulge_diskins, mbulge_mergers_assembly, mbulge_diskins_assembly,
     sAM_atomic_disk, vmax, rgas_disk) = hdf5_data

    bin_it = functools.partial(us.wmedians, xbins=xmf, low_numbers=False, nmin=50)

    ind = np.where(sfr_disk+sfr_burst <= 0)
    sfr_disk[ind] = 1e-10
    ind = np.where((sfr_disk+sfr_burst > 0) & (mdisk+mbulge > 0) & ((sfr_disk+sfr_burst)/(mdisk+mbulge) > 0))
    mainseqsf[index,:] = bin_it(x=np.log10((mdisk[ind]+mbulge[ind])/h0), y=np.log10((sfr_disk[ind]+sfr_burst[ind])/h0/GyrToYr))

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist = [0]
    #zlist = (0.005, 0.2, 0.5 , 0.8 , 1.1 , 1.5 , 2.2 , 2.9 , 3.9, 5.1)

    plt = common.load_matplotlib()

    mainseq     = np.zeros(shape = (len(zlist), 3, len(xmf)))


    fields = {'galaxies': ('sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge',
                           'mgas_metals_disk', 'mgas_metals_bulge',
                           'mstars_metals_disk', 'mstars_metals_bulge', 'type', 
                           'mvir_hosthalo', 'rstar_bulge', 'mstars_burst_mergers', 
                           'mstars_burst_diskinstabilities', 'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly',
                           'specific_angular_momentum_disk_gas_atom', 'vmax_subhalo', 'rgas_disk')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        prepare_data(hdf5_data, index, mainseq)

        h0 = hdf5_data[0]
        volh = hdf5_data[1]
        if index == 0:
            (sfr_disk, sfr_burst, mdisk, mbulge) = hdf5_data[2:6]
            sfr_seq = np.zeros(shape = (2, len(mdisk)))
            ind = np.where(sfr_disk + sfr_burst <= 0)
            sfr_disk[ind] = 2e-10 * GyrToYr * h0 #assume a minimum

            ind  = np.where((sfr_disk + sfr_burst > 0) & (mdisk + mbulge > 0))
            sfr_seq[0,ind] = np.log10((mdisk[ind]+mbulge[ind])/h0)
            sfr_seq[1,ind] = np.log10((sfr_disk[ind] + sfr_burst[ind]) / h0 / GyrToYr)

    plot_sfr_mstars_z0(plt, outdir, obsdir, h0, sfr_seq, mainseq)


if __name__ == '__main__':
    main(*common.parse_args())
