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

import numpy as np

import common
import utilities_statistics as us

mlow = 10
mupp = 15
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

def prepare_data(hdf5_data, index, massstellarh, massstellarh_v2):

    Omegab = 0.0491
    OmegaM = 0.3121

    (h0, _, mdisk, mbulge, ms_halo, mhalo, typeg, ms_tidally_stripped, id_halo) = hdf5_data


    ids_halo = np.unique(id_halo)

    mshalo_gals = np.zeros(shape = (3, len(ids_halo)))
    for i,g in enumerate(ids_halo): 
        ind = np.where(id_halo == g)
        mshalo_gals[0,i] = np.sum(ms_tidally_stripped[ind])
        mshalo_gals[1,i] = np.sum(mdisk[ind]+mbulge[ind])
        ind = np.where( (id_halo == g) & (typeg == 0))
        mshalo_gals[2,i] = mhalo[ind]

    ind = np.where((typeg <= 0) & (ms_halo > 0))
    massstellarh[index,:] = us.wmedians(x=np.log10(mhalo[ind]) - np.log10(float(h0)),
                                   y=np.log10(ms_halo[ind]) - np.log10(mhalo[ind]),
                                   xbins=xmf)
    ind = np.where(mshalo_gals[0,:] > 0)
    massstellarh_v2[index,:] = us.wmedians(x=np.log10(mshalo_gals[2,ind]) - np.log10(float(h0)),
                                   y=np.log10(mshalo_gals[0,ind]) - np.log10(mshalo_gals[1,ind]),
                                   xbins=xmf)

def plot_stellar_halo_z(plt, outdir, massstellarh, massstellarh_v2, zlist):

    fig = plt.figure(figsize=(7, 5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo, DM}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm M_{\\star,halo}/M_{\\star,total})$"
    xmin, xmax, ymin, ymax = 10.5, 15, -4, 0

    # z=0 ##################################
    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    plt.subplots_adjust(left=0.2)

    def compute_points(massstellarh):
        ind = np.where(massstellarh[i,0,:] != 0)
        xplot = xmf[ind]
        yplot = massstellarh[i,0,ind]
        errdn = massstellarh[i,1,ind]
        errup = massstellarh[i,2,ind]
        return (xplot, yplot, errdn, errup)
    
    cols = ['DarkRed', 'DarkOrange', 'LimeGreen', 'Turquoise', 'SteelBlue', 'Purple']
    #Predicted stellar-halo mass relation
    for i,z in enumerate(zlist):
        #(xplot, yplot, errdn, errup) = compute_points(massstellarh)
        #ax.plot(xplot, yplot[0], linestyle = 'solid', color=cols[i], label = 'z=%s' % str(z))
        #ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=cols[i], alpha=0.5, interpolate=True)
        #ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=cols[i], alpha=0.5, interpolate=True)

        (xplot, yplot, errdn, errup) = compute_points(massstellarh_v2)
        ax.plot(xplot, yplot[0], linestyle = 'solid', color=cols[i], label = 'z=%s' % str(z))
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=cols[i], alpha=0.5, interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=cols[i], alpha=0.5, interpolate=True)

    common.prepare_legend(ax, cols, loc=2)


    common.savefig(outdir, fig, 'stellar_halo_mass_z.pdf')


def main(modeldir, outdir, redshift_table, subvols):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstellar_halo', 'mvir_hosthalo',
                           'type','mstars_tidally_stripped','id_halo')}

    zlist = (0, 0.5, 1, 2, 3, 4)
    snapshots = redshift_table[zlist]
    massstellarh = np.zeros(shape = (len(zlist), 3, len(xmf)))
    massstellarh_v2 = np.zeros(shape = (len(zlist), 3, len(xmf)))

    for idx, snapshot in enumerate(snapshots):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        prepare_data(hdf5_data, idx, massstellarh, massstellarh_v2)

    plot_stellar_halo_z(plt, outdir, massstellarh, massstellarh_v2, zlist)


if __name__ == '__main__':
    main(*common.parse_args(requires_observations=False))
