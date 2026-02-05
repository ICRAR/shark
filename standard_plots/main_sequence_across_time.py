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
mlow = 7.75
mupp = 12
dm = 0.5
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

def plot_redshift_ms(ax,m,sfr,sfrs):
    sfr_hist_active = np.loadtxt("SFHs_ActiveGalaxies_PMill.txt")
    ms_hist_active = np.loadtxt("Masshistory_ActiveGalaxies_PMill.txt")
    lbts = np.loadtxt("LBT_PMill.txt")
    zs = us.redshift(lbts, h=0.6751, omegam=0.3121, omegal=0.6879)
    ngals = len(sfr_hist_active[:,0])
    nlbts = len(sfr_hist_active[0,:])
    #put histories in single vectors 
    sfr_all = np.zeros(shape = ngals * nlbts)
    ms_all = np.zeros(shape = ngals * nlbts)
    lbt_all = np.zeros(shape = ngals * nlbts)
  
    p = 0
    for g in range(0,ngals):
        for s in range(0,nlbts):
            sfr_all[p] = sfr_hist_active[g,s]
            ms_all[p] = ms_hist_active[g,s]
            lbt_all[p] = zs[s]
            p = p + 1
  
    binsz = np.array([15,9.5,8.5,7.5,6.5,5.5,4.5,3.5,2.5])
    cols = ['DarkSlateGray','DimGray','SlateGray','Gray','LightSlateGray','DarkGray','Silver','LightGray','Gainsboro']
    for i in range(len(m)-1):
        ind = np.where((ms_all >= 10**m[i]) & (ms_all < 10**m[i+1]))
        zplot = np.median(lbt_all[ind])
        print(zplot)
        select = np.where(binsz >= zplot)
        ncolsin = len(binsz[select])
        colin = cols[select[0][ncolsin-1]]
        ax.fill_between([m[i],m[i+1]],[sfr[i]-sfrs[i],sfr[i+1]-sfrs[i+1]],[sfr[i]+sfrs[i],sfr[i+1]+sfrs[i+1]],color=colin)
    #make the labels
    xs = 8.1
    dx = 0.25
    ys = 2.3
    dy = 0.15
    labels = [10, 9, 8, 7, 6, 5, 4, 3]
    for j in range(0,len(binsz)-1):
        ax.fill_between([xs + j * dx, xs + (j + 1) * dx], [ys, ys], [ys + dy, ys + dy], color=cols[j])
        ax.text(xs + j * dx + dx * 0.3, ys + dy + 0.02, str(labels[j]))
    ax.text(8.7, ys + dy + 0.15, 'median redshift')
  


def plot_sfr_mstars_evolution(plt, outdir, obsdir, zplots, main_seq):

    bin_it = functools.partial(us.wmedians, xbins=xmf, nmin=50)

    fig = plt.figure(figsize=(7,6))
    ytit = "$\\rm log_{10} (\\rm SFR/M_{\\odot})\\, yr^{-1})$"
    xtit = "$\\rm log_{10} (M_{\\star}/\\rm M_{\\odot})$"
    xmin, xmax, ymin, ymax = 8, 11.5, -0.4, 2.7
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    lws = [1,2, 3,5]
    ls = ['solid', 'solid', 'dashed', 'dotted']
    zplots = [9,7,5,4]

    for i in range(0,len(zplots)):
        #Predicted relation
        ind = np.where(main_seq[i,0,:] != 0)
        yplot = main_seq[i,0,ind]
        xplot = xmf[ind]
        errdn =  main_seq[i,1,ind]
        errup =  main_seq[i,2,ind]

        ax.plot(xplot,yplot[0],ls=ls[i], color='DarkOrange', lw=lws[i], label="z=%s" % str(zplots[i]))
        ax.fill_between(xplot,yplot[0]-errdn[0], yplot[0]+errup[0],ls=ls[i], color='DarkOrange', lw=lws[i], alpha=0.2)

    sm, sfrms, sfrs = np.loadtxt("MainSequenceActive.txt", unpack = True, usecols = [0,1,2])
    plot_redshift_ms(ax,sm,sfrms,sfrs)
    ax.plot(sm,sfrms,linestyle='solid', color='k', label = "control MS galaxies")
    ax.plot(sm,sfrms - sfrs,linestyle='dotted', color='k')
    ax.plot(sm,sfrms + sfrs,linestyle='dotted', color='k')

    common.prepare_legend(ax, ['DarkOrange','DarkOrange', 'DarkOrange', 'DarkOrange', 'k'], loc = 4)
    ax.text(8.2, 2,'SHARK', fontsize=16)

    # Legend
    plt.tight_layout()
    common.savefig(outdir, fig, 'main_sequence_evolution_shark.pdf')


def prepare_data(hdf5_data, index, mainseqsf):

    (h0, volh, sfr_disk, sfr_burst, mdisk, mbulge, typeg) = hdf5_data

    bin_it = functools.partial(us.wmedians, xbins=xmf, nmin=10)

    ssfr_in = (sfr_disk + sfr_burst)/(mdisk + mbulge)/1e9 

    ind = np.where((mdisk + mbulge > 1e8*h0) & (ssfr_in > 1e-9))
    mainseqsf[index,:] = bin_it(x=np.log10((mdisk[ind]+mbulge[ind])/h0), y=np.log10((sfr_disk[ind]+sfr_burst[ind])/h0/GyrToYr))

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist = [9,7,5,4]
    #zlist = (0.005, 0.2, 0.5 , 0.8 , 1.1 , 1.5 , 2.2 , 2.9 , 3.9, 5.1)

    plt = common.load_matplotlib()

    mainseq     = np.zeros(shape = (len(zlist), 3, len(xmf)))

    fields = {'galaxies': ('sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge','type')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        prepare_data(hdf5_data, index, mainseq)

    plot_sfr_mstars_evolution(plt, outdir, obsdir, zlist, mainseq)


if __name__ == '__main__':
    main(*common.parse_args())
