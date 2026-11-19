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

import functools
import numpy as np
import os

import common
import utilities_statistics as us

##################################

# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3
c_light = 299792458.0 #m/s

# Mass function initialization

mlow = -30 + 5.0 * np.log10(0.677)
mupp = -10 + 5.0 * np.log10(0.677)
dm = 0.5
mbins = np.arange(mlow,mupp,dm)
xlf   = mbins + dm/2.0


mflow = 5
mfupp = 14
dmf = 0.5
mfbins = np.arange(mflow,mfupp,dmf)
xmf = mfbins + dmf/2.0

mflow2 = 10
mfupp2 = 14
dmf2 = 0.25
mfbins2 = np.arange(mflow2,mfupp2,dmf2)
xmf2 = mfbins2 + dmf2/2.0


def plot_mass_evo(plt, outdir, obsdir, mstarbeta, z):

    fig = plt.figure(figsize=(8,4))
    xmin, xmax, ymin, ymax = 5, 11, 0.1, 3e3

    ytit="$\\rm \\dot{M}_{\\rm out}/SFR$"
    xtit="$\\rm log_{10}(M_{\\star}/M_{\\odot})$"

    ax = fig.add_subplot(111)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1e2, 1e2))
    ax.set_yscale('log')
    cols = ('k','Indigo','purple','Navy','Aquamarine', 'DarkOliveGreen', 'Green','Gold','red','DarkRed') 


    idxs = range(0,len(z))

    for idx in range(0,len(z)):
  
        #Predicted LF
        ind = np.where(mstarbeta[idx,0,:] != 0)
        print(ind)
        x = xmf[ind]
        y = mstarbeta[idx,0,ind]
        ydn = mstarbeta[idx,1,ind]
        yup = mstarbeta[idx,2,ind]
        ax.plot(x, y[0],color=cols[idx], linestyle='solid', linewidth = 4, label="z=%s" % str(z[idx]))
        ax.plot(x, y[0] - ydn[0],color=cols[idx], linestyle='dotted', linewidth = 4)
        ax.plot(x, y[0] + yup[0],color=cols[idx], linestyle='dotted', linewidth = 4)

    common.prepare_legend(ax, cols, bbox_to_anchor=(1.11, -0.025))

    plt.tight_layout()
    common.savefig(outdir, fig, "StellarMass_MassLoading.pdf")

def prepare_data(hdf5_data, index, mstarbeta, zlist):
  
    bin_it = functools.partial(us.wmedians, xbins=xmf, nmin=10)
    bin_it2 = functools.partial(us.wmedians, xbins=xmf2, nmin=10)

    #star_formation_histories and SharkSED have the same number of galaxies in the same order, and so we can safely assume that to be the case.
    #to select the same galaxies in galaxies.hdf5 we need to ask for all of those that have a stellar mass > 0, and then assume that they are in the same order.

    (h0, volh, mdisk, mbulge, mhalo, mshalo, typeg, age,
     sfr_disk, sfr_burst, id_gal, vvir_halo) = hdf5_data

    redshift_power = 0.13405459588358698
    v_sn = 120
    beta_disk = 3.79746174188
    eps_halo = 2.0
    min_beta = 0.104050197191
   
    age_univ =  us.look_back_time(zlist[index])
    vhot = v_sn * (age_univ)**redshift_power
    const_sn =  (vhot/vvir_halo)**beta_disk
    ind = np.where(const_sn < min_beta)
    const_sn[ind] = min_beta

    ind = np.where(mdisk + mbulge > 0)
    mstarbeta[index,:] = bin_it(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=const_sn[ind])

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()

    Variable_Ext = True
    multiple_batches = False
    #False

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age',
                           'sfr_disk', 'sfr_burst', 'id_galaxy', 'vvir_subhalo')}


    #z = (6, 6, 6.0, 8.0, 10.0, 13.0, 15.0)  #(8.0, 9, 10.5, 12.5) #, 1.0, 1.5, 2.0)
    z = (0.1, 1 , 2, 4, 6.0, 8.0, 10.0, 12.7, 15.0, 17.0)  #(8.0, 9, 10.5, 12.5) #, 1.0, 1.5, 2.0)
    #z = (17.0, 17.0, 17.0, 17.0, 17.0, 17.0, 17.0)  #(8.0, 9, 10.5, 12.5) #, 1.0, 1.5, 2.0)

    snapshots = redshift_table[z]
    mstarbeta = np.zeros(shape = (len(z), 3, len(xmf)))

    # Create histogram
    for index, snapshot in enumerate(snapshots):
        if(multiple_batches == False):
            hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        else:
            hdf5_data = common.read_data_multiple_batches(model_dir, snapshot, fields)

        prepare_data(hdf5_data, index, mstarbeta, z)

        h0, volh = hdf5_data[0], hdf5_data[1]
    plot_mass_evo(plt, outdir, obsdir, mstarbeta, z)

if __name__ == '__main__':
    main(*common.parse_args())
