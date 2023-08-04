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
from scipy import interpolate
import common
import os
import math 
import h5py

##################################

mlow = 5
mupp = 14
dm = 0.125
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

def plot_COLF_z2(plt, outdir, obsdir, hist_COLF):

    fig = plt.figure(figsize=(5,4.5))
    xtit = "$\\rm log_{10} (\\rm L^{\\prime}\\, [K\\, km/s\\, pc^2])$"
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}(\\rm L^{\\prime})/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 6.1, 12, -6.5, -0.8
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    z = [2, 2.5, 3]
    cols = ['red', 'blue', 'darkgreen']
    labels = ['CO(1-0)', 'CO(2-1)', 'CO(3-2)']
    linestyles = ['dotted', 'solid','dashed']
    for i,red in enumerate(z):
        for j in range(0,3):
           y = hist_COLF[i,j,:]
           ind = np.where(y < 0.)
           ax.plot(xmf[ind],y[ind], color=cols[j], linestyle=linestyles[i], label=labels[j] + " z=%s" % str(red))

    common.prepare_legend(ax, ['red', 'blue', 'darkgreen', 'red', 'blue', 'darkgreen', 'red', 'blue', 'darkgreen'])
    common.savefig(outdir, fig, 'COLF_z2p5.pdf')


def prepare_data(hdf5_data, hdf5_data_co, index, model_dir, snapshot, obsdir, hist_COLF):

    (h0, volh, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly,
     mBH, rdisk, rbulge, typeg, specific_angular_momentum_disk_star, specific_angular_momentum_bulge_star,
     specific_angular_momentum_disk_gas, specific_angular_momentum_bulge_gas, specific_angular_momentum_disk_gas_atom,
     specific_angular_momentum_disk_gas_mol, lambda_sub, mvir_s, mgas_disk, mgas_bulge, matom_disk, mmol_disk, matom_bulge,
     mmol_bulge, mbh_acc_hh, mbh_acc_sb, rgas_disk, rgas_bulge) = hdf5_data

    (lco_disk, lco_bulge) = hdf5_data_co
    co_total_10 = (lco_disk[:,0] + lco_bulge[:,0]) * 3.25e7 / (115.2712018)**2.0 / (4.0*np.pi) #in K km/s pc^2
    co_total_21 = (lco_disk[:,1] + lco_bulge[:,1]) * 3.25e7 / (230.5380000)**2.0 / (4.0*np.pi) #in K km/s pc^2
    co_total_32 = (lco_disk[:,2] + lco_bulge[:,2]) * 3.25e7 / (345.7959899)**2.0 / (4.0*np.pi) #in K km/s pc^2

    ind = np.where((np.log10(co_total_10) > 6) & (np.log10(co_total_10) < 15))
    H_COLF, _ = np.histogram(np.log10(co_total_10[ind]),bins=np.append(mbins,mupp))
    hist_COLF[index,0,:] = hist_COLF[index,0,:] + H_COLF

    ind = np.where((np.log10(co_total_21) > 6) & (np.log10(co_total_21) < 15))
    H_COLF, _ = np.histogram(np.log10(co_total_21[ind]),bins=np.append(mbins,mupp))
    hist_COLF[index,1,:] = hist_COLF[index,1,:] + H_COLF

    ind = np.where((np.log10(co_total_32) > 6) & (np.log10(co_total_32) < 15))
    H_COLF, _ = np.histogram(np.log10(co_total_32[ind]),bins=np.append(mbins,mupp))
    hist_COLF[index,2,:] = hist_COLF[index,2,:] + H_COLF

    return(h0, volh)

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):

    plt = common.load_matplotlib()

    zlist = (2, 2.5, 3)

    hist_COLF = np.zeros(shape = (len(zlist), 3, len(mbins)))

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'm_bh', 'rstar_disk', 'rstar_bulge', 'type',
                           'specific_angular_momentum_disk_star', 'specific_angular_momentum_bulge_star',
                           'specific_angular_momentum_disk_gas', 'specific_angular_momentum_bulge_gas',
                           'specific_angular_momentum_disk_gas_atom', 'specific_angular_momentum_disk_gas_mol',
                           'lambda_subhalo', 'mvir_subhalo', 'mgas_disk', 'mgas_bulge','matom_disk', 'mmol_disk',
                           'matom_bulge', 'mmol_bulge', 'bh_accretion_rate_hh', 'bh_accretion_rate_sb','rgas_disk', 'rgas_bulge')}
    fields_co = {'galaxies': ('LCO_disk','LCO_bulge')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        hdf5_data_co = common.read_co_data(model_dir, snapshot, fields_co, subvols)

        (h0, volh) = prepare_data(hdf5_data, hdf5_data_co, index, model_dir, snapshot, obs_dir, hist_COLF)

    ind = np.where(hist_COLF[:] != 0)
    hist_COLF[ind] = np.log10(hist_COLF[ind]/volh * h0**3/ dm)
 
    for i,z in enumerate(zlist):
        print("#CO LFs at redshift", z)
        print("log10(L'CO[K km/s pc^2]) log10(phi(1-0)[Mpc^-3 dex^-1]) log10(phi(2-1)[Mpc^-3 dex^-1]) log10(phi(3-2)[Mpc^-3 dex^-1])")
        for a,b,c,d in zip(xmf, hist_COLF[i,0,:], hist_COLF[i,1,:], hist_COLF[i,2,:]):
            print(a,b,c,d)

    plot_COLF_z2(plt, output_dir, obs_dir, hist_COLF)

if __name__ == '__main__':
    main(*common.parse_args())
