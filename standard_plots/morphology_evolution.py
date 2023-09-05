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
zlist_old = xz

zlist= np.zeros(shape = len(zlist_old) + 1)
zlist[0] = 0
zlist[1:len(zlist)] = zlist_old[:]

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


def prepare_data(hdf5_data, index, BT_fractions, BT_fractions_distribution, mass_by_comp):

    (h0, volh, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     mBH, rdisk, rbulge, typeg, specific_angular_momentum_disk_star, specific_angular_momentum_bulge_star, 
     specific_angular_momentum_disk_gas, specific_angular_momentum_bulge_gas, specific_angular_momentum_disk_gas_atom, 
     specific_angular_momentum_disk_gas_mol, lambda_sub, mvir_s, mgas_disk, mgas_bulge, matom_disk, mmol_disk, matom_bulge, 
     mmol_bulge, mbh_acc_hh, mbh_acc_sb) = hdf5_data

    volume = volh / h0**3
    print("volume", volume)
    mstars_tot = (mdisk+mbulge)/h0
    mstars_disk_tot = np.sum(mdisk/h0) / volume
    mstars_bulge_mergers_tot = np.sum(mburst_mergers + mstars_bulge_mergers_assembly)/h0/volume
    mstars_bulge_diskins_tot =np.sum(mburst_diskins + mstars_bulge_diskins_assembly)/h0/volume

    mass_by_comp[index,0] = mstars_disk_tot
    mass_by_comp[index,1] = mstars_bulge_mergers_tot
    mass_by_comp[index,2] = mstars_bulge_diskins_tot

    bt = mbulge / (mdisk+mbulge)

    mall = np.sum(mstars_tot)
    BT_fractions[index, 0] = np.sum(mall)/volh * h0**2.0
    for i in range(0,len(btbins)-1):
        ind = np.where( ( mstars_tot > 0) & (bt >= btbins[i]) & (bt <= btbins[i+1]))
        BT_fractions[index, i+1] = np.sum(mstars_tot[ind])/volh * h0**2.0

    ind = np.where(mstars_tot/h0 > 1e9)
    H, _ = np.histogram(bt[ind],bins=np.append(btbins2, btupp))
    BT_fractions_distribution[0,index,:] = BT_fractions_distribution[0,index,:] + H
    #normalize to get PDF
    BT_fractions_distribution[0,index,:] = BT_fractions_distribution[0,index,:] / (len(bt[ind]) * dbt)

def plot_bt_fractions(plt, outdir, obsdir, BT_fractions, BT_fractions_distribution):

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

    BT_fractions = np.zeros(shape = (len(zlist), len(btbins)))
    BT_fractions_distribution = np.zeros(shape = (2, len(zlist), len(btbins2)))
    mass_by_comp = np.zeros(shape = (len(zlist), 3))

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        prepare_data(hdf5_data, index, BT_fractions, BT_fractions_distribution, mass_by_comp)

    plot_bt_fractions(plt, outdir, obsdir, BT_fractions, BT_fractions_distribution)

    for a,b,c,d,e in zip(zlist, us.look_back_time(zlist), mass_by_comp[:,0], mass_by_comp[:,1],  mass_by_comp[:,2]):
        print(a,b,c,d, e)

if __name__ == '__main__':
    main(*common.parse_args())
