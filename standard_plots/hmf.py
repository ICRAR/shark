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
import utilities_statistics as us

##################################
# Mass function initialization
mlow = 10
mupp = 15
dm = 0.3
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0 

smlow = 5
smupp = 14
dsm   = 0.2
smbins = np.arange(smlow,smupp,dsm)
xsmf   = smbins + dsm/2.0

Omegab = 0.0491
OmegaM = 0.3121
fb     = 0.0491/(0.3121-0.0491)

def plot_halomf_z(plt, outdir, obsdir, z, h0, hist, histsh, plotz):

    xtit="$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit="$\\rm log_{10}(\Phi/dlog{\\rm M_{\\rm halo}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 10.1, 15, -6, -1
    xleg = xmax - 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    fig = plt.figure(figsize=(7,7))

    subplots = (321, 322, 323, 324, 325, 326)
    idx = (0, 1, 2, 3, 4, 5)
    for subplot, idx, z, plot_this_z in zip(subplots, idx, z, plotz):

        ax = fig.add_subplot(subplot)
        if (idx == 0 or idx == 2 or idx == 4):
            ytitplot = ytit
        else:
            ytitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitplot, locators=(0.1, 1, 0.1))
        ax.text(xleg,yleg, 'z=%s' % (str(z)))

        if(idx < 4):
           #HMF calc HMF calculated by Sheth & Tormen (2001)
           lmp, dp = common.load_observation(obsdir, 'mf/HMF/mVector_PLANCK-SMT_z%s.dat' % str(z).replace('.', ''), [0, 7])
           lmp_plot = np.log10(lmp) - np.log10(h0)
           dp_plot = np.log10(dp) + np.log10(pow(h0,3.))
           if idx == 0:
               ax.plot(lmp_plot,dp_plot,'b', label = 'HMF calc')
           elif idx > 0:
               ax.plot(lmp_plot,dp_plot,'b')
   
        #Predicted HMF
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

    xtit="$\\rm log_{10}(\\rm M/M_{\odot})$"
    ytit="$\\rm log_{10}(\Phi/dlog{\\rm M}/{\\rm Mpc}^{-3})$"
    xmin, xmax, ymin, ymax = 8, 15, -6, 1.2
    xleg = xmax - 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    fig = plt.figure(figsize=(5,5))

    idx = 0
    ax = fig.add_subplot(111)
    ytitplot = ytit
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitplot, locators=(0.1, 1, 0.1))

    #lmp, dp = common.load_observation(obsdir, 'mf/HMF/mVector_PLANCK-SMT_z0_extended.dat', [0, 7])
    #lmp_plot = np.log10(lmp) - np.log10(h0)
    #dp_plot = np.log10(dp) + np.log10(pow(h0,3.))
    #ax.plot(lmp_plot,dp_plot,'k')
    #ax.plot(lmp_plot+np.log10(fb),dp_plot,'b', linestyle='dashed')

    y = hist[idx,:]
    ind = np.where((y < 0.) & (xmf > 10.3))
    ax.plot(xmf[ind],y[ind],'k')
    ax.plot(xmf[ind]+np.log10(fb),y[ind],'b', linestyle='dashed')

    lm, p, dpdn, dpup = common.load_observation(obsdir, 'mf/SMF/GAMAII_BBD_GSMFs.dat', [0,1,2,3])
    xobs = lm
    indx = np.where(p > 0)
    yobs = np.log10(p[indx])
    ydn = yobs - np.log10(p[indx]-dpdn[indx])
    yup = np.log10(p[indx]+dpup[indx]) - yobs
    ax.errorbar(xobs[indx], yobs, ydn, yup, 'ro', label='Wright+17')


    smfdensity = common.load_observation(obsdir, 'Models/SharkVariations/SMF_FeedbackExperiment.dat', [0])
    ynofeed = smfdensity[0:len(xsmf)-1]
    yreio   = smfdensity[len(xsmf):2*len(xsmf)-1]
    ystarf  = smfdensity[2*len(xsmf):3*len(xsmf)-1]
    yfinal  = smfdensity[3*len(xsmf):4*len(xsmf)-1]
    
    ind = np.where(ynofeed != 0)
    ax.plot(xsmf[ind],ynofeed[ind],'DarkRed')
    ind = np.where(ystarf != 0)
    ax.plot(xsmf[ind],ystarf[ind],'LightCoral')
    ind = np.where(yfinal != 0)
    ax.plot(xsmf[ind],yfinal[ind],'Orange')

    common.prepare_legend(ax, ['r'])

    common.savefig(outdir, fig, "halomf_z0.pdf")


def plot_concentration_z(plt, outdir, z, concentration):

    xtit = "$\\rm log_{10}(M_{\\rm halo}/M_{\\odot})$"
    ytit = "$\\rm log_{10}(c_{\\rm vir})$"

    xmin, xmax, ymin, ymax = 10.1, 15, 0, 1.5
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    fig = plt.figure(figsize=(7, 7))

    subplots = (321, 322, 323, 324, 325, 326)
    idx = (0, 1, 2, 3, 4, 5)

    for subplot, idx, redshift in zip(subplots, idx, z):

        ax = fig.add_subplot(subplot)

        if (idx == 0 or idx == 2 or idx == 4):
            ytitplot = ytit
        else:
            ytitplot = ' '

        common.prepare_ax(
            ax, xmin, xmax, ymin, ymax,
            xtit, ytitplot,
            locators=(1, 1, 0.25, 0.25)
        )

        ax.text(xleg, yleg, 'z=%s' % str(redshift))

        # Median concentration and 16th-84th percentile scatter
        y = concentration[idx, 0, :]
        ydn = concentration[idx, 1, :]
        yup = concentration[idx, 2, :]

        ind = np.where(y > 0)

        ax.plot(xmf[ind], y[ind], 'r')

        ax.fill_between(
            xmf[ind],
            y[ind] - ydn[ind],
            y[ind] + yup[ind],
            color='r',
            alpha=0.2
        )

    common.savefig(outdir, fig, "halo_concentration_z.pdf")

def prepare_data(hdf5_data, hist, histsh, concentration, index):

    (h0, _, mdisk, mbulge, mhalo, mshalo, typeg, 
    mhalo_all, cnfw) = hdf5_data

    mass   = np.zeros(shape = len(mhalo))
    masssh = np.zeros(shape = len(mhalo))

    # halo mass function calculation
    ind = np.where((typeg <= 0) & (mdisk+mbulge > 1e5))
    mass[ind] = np.log10(mhalo[ind]) - np.log10(float(h0))
    ind = np.where((typeg <= 1) & (mdisk+mbulge > 1e5))
    masssh[ind] = np.log10(mshalo[ind]) - np.log10(float(h0))

    H, bins_edges = np.histogram(mass,bins=np.append(mbins,mupp))
    hist[index,:] = hist[index,:] + H
    H, bins_edges = np.histogram(masssh,bins=np.append(mbins,mupp))
    histsh[index,:] = histsh[index,:] + H

    # Halo concentration-mass relation
    ind = np.where(
        (mhalo_all > 0) &
        (cnfw > 0) &
        np.isfinite(mhalo_all) &
        np.isfinite(cnfw)
    )

    logm = np.log10(mhalo_all[ind]) - np.log10(float(h0))
    logc = np.log10(cnfw[ind])

    concentration[index, :, :] = us.wmedians(x=logm, y=logc,xbins=xmf)

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type'),
                'halo': ('mvir','concentration')}

    z = (0, 0.5, 1, 2, 4, 6)
    snapshots = redshift_table[z]

    # Create histogram
    hist = np.zeros(shape = (len(z), len(mbins)))
    histsh = np.zeros(shape = (len(z), len(mbins)))
    plotz = np.empty(shape=(len(z)), dtype=np.bool_)
    concentration = np.zeros(shape=(len(z), 3, len(mbins)))

    for index, snapshot in enumerate(snapshots):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        prepare_data(hdf5_data, hist, histsh, concentration, index)

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

    plot_halomf_z(plt, outdir, obsdir, z, h0, hist, histsh, plotz)
    plot_concentration_z(plt, outdir, z, concentration)

if __name__ == '__main__':
    main(*common.parse_args())
