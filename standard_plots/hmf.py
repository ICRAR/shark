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

import numpy as np

import common


##################################
# Mass function initialization
mlow = 10
mupp = 15
dm = 0.3
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0 

def plot_halomf_z(plt, outdir, obsdir, h0, hist, histsh, plotz):

    xtit="$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit="$\\rm log_{10}(\Phi/dlog{\\rm M_{\\rm halo}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 10.1, 15, -6, -1
    xleg = xmax - 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    fig = plt.figure(figsize=(7,7))

    subplots = (221, 222, 223, 224)
    idx = (0, 1, 2, 3)
    z = (0, 0.5, 1, 2)
    for subplot, idx, z, plot_this_z in zip(subplots, idx, z, plotz):

        ax = fig.add_subplot(subplot)
        if (idx == 0 or idx == 2):
            ytitplot = ytit
        else:
            ytitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitplot, locators=(0.1, 1, 0.1))
        ax.text(xleg,yleg, 'z=%s' % (str(z)))

        #HMF calc HMF calculated by Sheth & Tormen (2001)
        lmp, dp = common.load_observation(obsdir, 'mf/HMF/mVector_PLANCK-SMT_z%s.dat' % str(z).replace('.', ''), [0, 7])
        lmp_plot = np.log10(lmp) - np.log10(h0)
        dp_plot = np.log10(dp) + np.log10(pow(h0,3.))
        if idx == 0:
            ax.plot(lmp_plot,dp_plot,'b', label = 'HMF calc')
        elif idx > 0:
            ax.plot(lmp_plot,dp_plot,'b')

        #Predicted HMF
        if plot_this_z:
            y = hist[idx,:]
            ind = np.where(y < 0.)
            if idx == 0:
                ax.plot(xmf[ind],y[ind],'r', label ='HMF Shark')
            if idx > 0:
                ax.plot(xmf[ind],y[ind],'r')
            y = histsh[idx,:]
            ind = np.where(y < 0.)
            if idx == 0:
                ax.plot(xmf[ind],y[ind],'r', linestyle='dashed', label ='SHMF Shark')
            if idx > 0:
                ax.plot(xmf[ind],y[ind],'r', linestyle='dashed')

        if idx == 0:
            common.prepare_legend(ax, ['b','r','r'])

    common.savefig(outdir, fig, "halomf_z.pdf")

def prepare_data(hdf5_data, hist, histsh, index):

    h0, _, mdisk, mbulge, mhalo, mshalo, typeg = hdf5_data
    mass   = np.zeros(shape = len(mhalo))
    masssh = np.zeros(shape = len(mhalo))

    ind = np.where((typeg <= 0) & (mdisk+mbulge > 1e5))
    mass[ind] = np.log10(mhalo[ind]) - np.log10(float(h0))
    ind = np.where((typeg <= 1) & (mdisk+mbulge > 1e5))
    masssh[ind] = np.log10(mshalo[ind]) - np.log10(float(h0))

    H, bins_edges = np.histogram(mass,bins=np.append(mbins,mupp))
    hist[index,:] = hist[index,:] + H
    H, bins_edges = np.histogram(masssh,bins=np.append(mbins,mupp))
    histsh[index,:] = histsh[index,:] + H

def main(model_dir, outdir, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type')}

    # Create histogram
    zlist = ["199","174", "156", "131"]
    hist = np.zeros(shape = (len(zlist), len(mbins)))
    histsh = np.zeros(shape = (len(zlist), len(mbins)))
    plotz = np.empty(shape=(len(zlist)), dtype=np.bool_)

    for index in range(0,4):

        hdf5_data = common.read_data(model_dir, zlist[index], fields, subvols)
        prepare_data(hdf5_data, hist, histsh, index)

        h0, volh = hdf5_data[0], hdf5_data[1]
        if(volh > 0.):
            vol = volh/pow(h0,3.)  # In Mpc^3
            hist[index,:]   = hist[index,:]/vol/dm
            histsh[index,:] = histsh[index,:]/vol/dm
            plotz[index]    = True
        else:
            plotz[index]    = False

    # Take logs
    ind = np.where(hist > 0.)
    hist[ind] = np.log10(hist[ind])

    ind = np.where(histsh > 0.)
    histsh[ind] = np.log10(histsh[ind])

    plot_halomf_z(plt, outdir, obsdir, h0, hist, histsh, plotz)

if __name__ == '__main__':
    main(*common.parse_args(requires_snapshot=False))
