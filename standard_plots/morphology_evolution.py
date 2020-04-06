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

# Initialize arguments
zlow = 0
zupp = 6
dz = 0.05
zbins = np.arange(zlow,zupp,dz)
xz   = zbins + dz/2.0
zlist = xz

##################################
#Constants
RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun

mlow = 6.5
mupp = 12.5
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

dmobs = 0.4
mbins_obs = np.arange(mlow,mupp,dmobs)
xmf_obs = mbins_obs + dmobs/2.0

vlow = 1.0
vupp = 3.0
dv   = 0.1
vbins = np.arange(vlow,vupp,dv)
xv    = vbins + dv/2.0

btbins = [0, 0.05, 0.95, 1]


def prepare_data(hdf5_data, index, BT_fractions):

    (h0, volh, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     mBH, rdisk, rbulge, typeg, specific_angular_momentum_disk_star, specific_angular_momentum_bulge_star, 
     specific_angular_momentum_disk_gas, specific_angular_momentum_bulge_gas, specific_angular_momentum_disk_gas_atom, 
     specific_angular_momentum_disk_gas_mol, lambda_sub, mvir_s, mgas_disk, mgas_bulge, matom_disk, mmol_disk, matom_bulge, 
     mmol_bulge, mbh_acc_hh, mbh_acc_sb) = hdf5_data

    mstars_tot = (mdisk+mbulge)/h0
    bt = mbulge / (mdisk+mbulge)
    ind = np.where(mstars_tot > 0)
    mall = np.sum(mstars_tot[ind])
    BT_fractions[index, 0] = np.sum(mall)/volh * h0**2.0
    for i in range(0,len(btbins)-1):
        ind = np.where( ( mstars_tot > 0) & (bt >= btbins[i]) & (bt <= btbins[i+1]))
        BT_fractions[index, i+1] = np.sum(mstars_tot[ind])/volh * h0**2.0

def plot_bt_fractions(plt, outdir, obsdir, BT_fractions):

    fig = plt.figure(figsize=(5,4.5))
    ytit = "$\\rm log_{10} (\\rm \\rho_{\\star}/ M_{\\odot} cMpc^{-3})$"
    xtit = "lookback time/Gyr"
    xmin, xmax, ymin, ymax = 0, 12, 6, 9
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #Predicted relation
    cols = ['k','DarkSlateBlue','GreenYellow', 'Crimson']
    labels = ['total', '$\\rm  B/T <0.05$', '$\\rm  0.05<B/T<0.95$', '$\\rm  B/T>0.95$']

    print (us.look_back_time(zlist), zlist)
    for i in range(0,len(cols)):
        ax.plot(us.look_back_time(zlist), np.log10(BT_fractions[:,i]), linestyle='solid',color=cols[i], label=labels[i])
        for a,b in zip(us.look_back_time(zlist), np.log10(BT_fractions[:,i])):
            print (a,b)
    common.prepare_legend(ax, cols, loc=4)
    common.savefig(outdir, fig, 'morphology_evolution.pdf')

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'm_bh', 'rstar_disk', 'rstar_bulge', 'type', 
                           'specific_angular_momentum_disk_star', 'specific_angular_momentum_bulge_star',
                           'specific_angular_momentum_disk_gas', 'specific_angular_momentum_bulge_gas',
                           'specific_angular_momentum_disk_gas_atom', 'specific_angular_momentum_disk_gas_mol',
                           'lambda_subhalo', 'mvir_subhalo', 'mgas_disk', 'mgas_bulge','matom_disk', 'mmol_disk', 
                           'matom_bulge', 'mmol_bulge', 'bh_accretion_rate_hh', 'bh_accretion_rate_sb')}

    BT_fractions = np.zeros(shape = (len(zlist), len(btbins)))
   
    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        prepare_data(hdf5_data, index, BT_fractions)

    plot_bt_fractions(plt, outdir, obsdir, BT_fractions)


if __name__ == '__main__':
    main(*common.parse_args())
