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

##################################
#Constants
RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun


slow = -12.0
supp = -8.0
ds = 0.2
sbins = np.arange(slow,supp,ds)
xsf = sbins + ds/2.0

def prepare_data(hdf5_data, index, hist_ssfr, mbins, dm):

    (h0, volh, mdisk, mbulge, sfr_disk, sfr_burst) = hdf5_data

    mstars_tot = (mdisk+mbulge)/h0
    sfr_tot = (sfr_disk + sfr_burst)/h0/1e9 #Msun/yr
    
    ind = np.where(mstars_tot > 0)
    mstars_tot[ind] = np.log10(mstars_tot[ind])
    ind = np.where(sfr_tot > 0)
    sfr_tot[ind] = np.log10(sfr_tot[ind])

    for i,m in enumerate(mbins):
        ind = np.where((mstars_tot > m-dm*0.5) & (mstars_tot<= m+dm*0.5))
        ssfr_in = sfr_tot[ind] - mstars_tot[ind]
        H, _ = np.histogram(ssfr_in,bins=np.append(sbins,supp))
        hist_ssfr[index,i,:] = hist_ssfr[index,i,:] + H

    vol = volh / h0**3
    return(vol)

def plot_ssfr(plt, outdir, obsdir, hist_ssfr):
    
    plots = [221, 222, 223, 224]

    fig = plt.figure(figsize=(10,10))
    ytit = "$\\rm log_{10} (\\rm \\rho_{\\rm sSFR}/ cMpc^{-3} dex^{-1})$"
    xtit = "$\\rm log_{10} (\\rm sSFR/yr^{-1})$"
    xmin, xmax, ymin, ymax = -12, -8, 1e-6, 0.1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    def plot_observations(ax,i):
        #load observations
        sm, phi1, phi1err, phi2, phi2err, phi3, phi3err, phi4, phi4err = common.load_observation(obsdir, 'SFR/SSFR_distributions_Katsianis.dat', [0,1,2,3,4,5,6,7,8])
        label = 'Katsianis+21'
        if i == 0:
           phi = phi1
           err = phi1err
        elif i == 1:
           phi = phi2
           err = phi2err
        elif i == 2:
           phi = phi3
           err = phi3err
        elif i == 3:
           phi = phi4
           err = phi4err
        ind = np.where(phi != 0)
        xplot = sm[ind]
        yplot = phi[ind]
        yerr = err[ind]
        if i == 0:
            ax.errorbar(xplot, yplot, yerr=yerr, ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D', label = label)
        else:
            ax.errorbar(xplot, yplot, yerr=yerr, ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D')

    labels = ['[9.5-10)','[10-10.5)','[10.5-11)','[11-11.5)']
    for i, s in enumerate(plots):
        ax = fig.add_subplot(s)
        if (i==1 or i == 3):
            ytitle = ' '
        else:
            ytitle = ytit
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle, locators=(0.1, 1, 0.1, 1))
        ax.set_yscale('log')
        z=0
        ind = np.where(hist_ssfr[z,i,:] != 0)
        xplot = xsf[ind]
        yplot = hist_ssfr[z,i,ind]
        ax.plot(xplot,yplot[0],'k', linestyle='solid', label="$\\rm log_{10}(M_{\\star}/M_{\\odot})=$" + labels[i])
        ax.plot(xplot-0.3,yplot[0],'k', linestyle='dashed', label="-0.3dex")

        plot_observations(ax,i)
        common.prepare_legend(ax, ['k','k'], loc=2)

    common.savefig(outdir, fig, 'ssfr_distribution_z0.pdf')

def main(modeldir, outdir, redshift_table, subvols, obsdir):


    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'sfr_disk', 'sfr_burst')}

    zlist = [0]
    mbins = [9.75,10.25,10.75,11.25]
    dm = 0.5

    hist_ssfr = np.zeros(shape = (len(zlist), len(mbins), len(sbins)))

    # Read data from each subvolume at a time and add it up
    # rather than appending it all together
    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        vol = prepare_data(hdf5_data, index, hist_ssfr, mbins, dm)

    hist_ssfr = hist_ssfr/vol/ds
    plot_ssfr(plt, outdir, obsdir, hist_ssfr)

if __name__ == '__main__':
    main(*common.parse_args())
