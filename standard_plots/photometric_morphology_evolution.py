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
zupp = 1
dz = 0.1
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

btlow = 0.0
btupp = 1.0
dbt   = 0.05
btbins2 = np.arange(btlow,btupp,dbt)
xbt    = btbins2 + dbt/2.0


def prepare_data(hdf5_data, seds, index, BT_fractions, BT_fractions_distribution):

    (h0, volh, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     mBH, rdisk, rbulge, typeg, specific_angular_momentum_disk_star, specific_angular_momentum_bulge_star, 
     specific_angular_momentum_disk_gas, specific_angular_momentum_bulge_gas, specific_angular_momentum_disk_gas_atom, 
     specific_angular_momentum_disk_gas_mol, lambda_sub, mvir_s, mgas_disk, mgas_bulge, matom_disk, mmol_disk, matom_bulge, 
     mmol_bulge, mbh_acc_hh, mbh_acc_sb) = hdf5_data

    
    mstars_tot = (mdisk+mbulge)/h0
    bt = mbulge / (mdisk+mbulge)
    ind = np.where(mstars_tot > 0)
    mstars = mstars_tot[ind]
    SEDs_dust_total = seds[1]
    SEDs_dust_bulge = seds[0]


    mall = np.sum(mstars_tot[ind])
    BT_fractions[index, 0] = np.sum(mall)/volh * h0**2.0
    for i in range(0,len(btbins)-1):
        ind = np.where( ( mstars_tot > 0) & (bt >= btbins[i]) & (bt <= btbins[i+1]))
        BT_fractions[index, i+1] = np.sum(mstars_tot[ind])/volh * h0**2.0

    ind = np.where(mstars_tot/h0 > 1e9)
    H, _ = np.histogram(bt[ind],bins=np.append(btbins2, btupp))
    BT_fractions_distribution[0,index,:] = BT_fractions_distribution[0,index,:] + H
    #normalize to get PDF
    BT_fractions_distribution[0,index,:] = BT_fractions_distribution[0,index,:] / (len(bt[ind]) * dbt)

    btr = 10.0**(SEDs_dust_bulge[4,:]/(-2.5)) / 10.0**(SEDs_dust_total[4,:]/(-2.5))
    ind = np.where(SEDs_dust_bulge[4,:] == -999)
    btr[ind] = 0
    ind = np.where(mstars > 1e9)
    H, _ = np.histogram(btr[ind],bins=np.append(btbins2, btupp))
    BT_fractions_distribution[1,index,:] = BT_fractions_distribution[1,index,:] + H
    #normalize to get PDF
    BT_fractions_distribution[1,index,:] = BT_fractions_distribution[1,index,:] / (len(mstars[ind]) * dbt)

def plot_bt_fractions(plt, outdir, obsdir, BT_fractions, BT_fractions_distribution):

    p= [0,1]
    names = ['stellar','rband']
    ytit = "PDF"
    xtit = "$B/T$"
    xmin, xmax, ymin, ymax = 0, 1, 0, 4
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    for j in p:
    
        fig = plt.figure(figsize=(5,4.5))
    
        ax = fig.add_subplot(111)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
        ax.text(xleg,yleg,names[j],fontsize=12)
 
        ind = np.where(zlist < 1.1)
        zlistin = zlist[ind]
        BT_fractions_distribution_ins = BT_fractions_distribution[j,ind,:]
        colors = plt.cm.jet(np.linspace(0,1,len(zlistin))) 
        #Predicted relation
        for i in range(0,len(zlistin)):
            ind = np.where(BT_fractions_distribution_ins[0,i,:] > 0)
            y = BT_fractions_distribution_ins[0,i,ind]
            if((i+1)%2 == 1):
               ax.plot(xbt[ind], y[0], linestyle='solid',color=colors[i], label='z=%s' % str(zlistin[i]))
            else:
               ax.plot(xbt[ind], y[0], linestyle='solid',color=colors[i])
    
        cols_for_label = colors[::2]
        common.prepare_legend(ax, cols_for_label, loc='upper center')
        common.savefig(outdir, fig, 'distribution_BT_evolution_'+names[j]+'.pdf')

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'm_bh', 'rstar_disk', 'rstar_bulge', 'type', 
                           'specific_angular_momentum_disk_star', 'specific_angular_momentum_bulge_star',
                           'specific_angular_momentum_disk_gas', 'specific_angular_momentum_bulge_gas',
                           'specific_angular_momentum_disk_gas_atom', 'specific_angular_momentum_disk_gas_mol',
                           'lambda_subhalo', 'mvir_subhalo', 'mgas_disk', 'mgas_bulge','matom_disk', 'mmol_disk', 
                           'matom_bulge', 'mmol_bulge', 'bh_accretion_rate_hh', 'bh_accretion_rate_sb')}

    fields_sed = {'SED/ab_dust': ('bulge_t','total'),}
    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5"

    BT_fractions = np.zeros(shape = (len(zlist), len(btbins)))
    BT_fractions_distribution = np.zeros(shape = (2, len(zlist), len(btbins2)))
  
    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        seds = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed, subvols, file_hdf5_sed)
        prepare_data(hdf5_data, seds, index, BT_fractions, BT_fractions_distribution)

    plot_bt_fractions(plt, outdir, obsdir, BT_fractions, BT_fractions_distribution)


if __name__ == '__main__':
    main(*common.parse_args())
