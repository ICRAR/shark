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

binsm = [8.5, 9.5, 10.5, 12.5]

recycle = 0.4588

def prepare_data(hdf5_data, sfh, hdf5_data_hist, LBT, delta_t):

    (h0, volh,  _, mstar, sfrdisk, sfrburst) = hdf5_data_hist
    (h0, volh, mdisk, mbulge) = hdf5_data
    (bulge_diskins_hist, bulge_mergers_hist, disk_hist) = sfh

    vol_sim = volh / h0**3
    print(sum(mdisk)/ sum(mdisk + mbulge))

    ms_true = mstar / h0 / (volh / h0**3.0)
    sfrall = (sfrdisk + sfrburst)/h0/ 1e9/(volh / h0**3.0)
    mstars_tot = (mdisk+mbulge)/h0
    bt = mbulge / (mdisk+mbulge)
    ind = np.where(mstars_tot > 0)
    mstars = np.log10(mstars_tot[ind])

    sfhs_masses = np.zeros(shape= (len(LBT), len(mbins)))
    sfhs_comp = np.zeros(shape= (len(LBT), 3))
    sfr_tot = np.zeros(shape= (len(LBT)))

    mshs_masses     = np.zeros(shape= (len(LBT), len(mbins)))
    mshs_comp = np.zeros(shape= (len(LBT), 3))
    mshs_masses_cum = np.zeros(shape= (len(LBT), len(mbins)))
    ms_tot = np.zeros(shape= (len(LBT)))
    ms_tot_bycomp = np.zeros(shape= (len(LBT), 3))
    ms_tot_bycomp_massivegals = np.zeros(shape= (len(LBT), 3))

    for i in range(0,len(binsm)):
        if ( i == 0 ):
             ind = np.where((mstars <= binsm[i]) & (mstars >= 7))
        else:
             ind = np.where((mstars <= binsm[i]) & (mstars > binsm[i-1]))
        for j in range(0,len(LBT)):
              sfhs_masses[j,i] = np.sum(disk_hist[ind,j] + bulge_mergers_hist[ind,j] + bulge_diskins_hist[ind,j]) / h0/ (volh / h0**3.0)
              mshs_masses[j,i] = np.sum(disk_hist[ind,j] + bulge_mergers_hist[ind,j] + bulge_diskins_hist[ind,j]) * delta_t[j] * 1e9 / h0 * (1.0 - recycle)

    for j in range(0,len(LBT)):
        sfr_tot[j] = np.sum(disk_hist[:,j] + bulge_mergers_hist[:,j] + bulge_diskins_hist[:,j]) / h0/ (volh / h0**3.0)
        for i in range(0,len(binsm)):
            mshs_masses_cum[j,i] = np.sum(mshs_masses[0:j,i]) / (volh / h0**3.0)
        sfhs_comp[j,0] = np.sum(disk_hist[:,j]) / h0/ (volh / h0**3.0)
        sfhs_comp[j,1] = np.sum(bulge_mergers_hist[:,j]) / h0/ (volh / h0**3.0)
        sfhs_comp[j,2] = np.sum(bulge_diskins_hist[:,j]) / h0/ (volh / h0**3.0)
        
        if j == 0:
           mshs_comp[j,0] = np.sum(disk_hist[:,j]) / h0/ (volh / h0**3.0) * 1e9 * delta_t[j] * (1.0 - recycle)
           mshs_comp[j,1] = np.sum(bulge_mergers_hist[:,j]) / h0/ (volh / h0**3.0) * 1e9  * delta_t[j] * (1.0 - recycle)
           mshs_comp[j,2] = np.sum(bulge_diskins_hist[:,j]) / h0/ (volh / h0**3.0) * 1e9  * delta_t[j] * (1.0 - recycle)
        else:
           mshs_comp[j,0] = np.sum(disk_hist[:,j]) / h0/ (volh / h0**3.0) * 1e9  * delta_t[j] * (1.0 - recycle) + mshs_comp[j-1,0]
           mshs_comp[j,1] = np.sum(bulge_mergers_hist[:,j]) / h0/ (volh / h0**3.0) * 1e9  * delta_t[j] * (1.0 - recycle) + mshs_comp[j-1,1]
           mshs_comp[j,2] = np.sum(bulge_diskins_hist[:,j]) / h0/ (volh / h0**3.0) * 1e9  * delta_t[j] * (1.0 - recycle) + mshs_comp[j-1,2] 

        ms_tot[j] = np.sum(mshs_masses_cum[j,:])
        ms_tot_bycomp[j,0] = np.sum(disk_hist[:,j]) / np.sum(disk_hist[:,j] + bulge_mergers_hist[:,j] + bulge_diskins_hist[:,j])
        ms_tot_bycomp[j,1] = np.sum(bulge_mergers_hist[:,j]) / np.sum(disk_hist[:,j] + bulge_mergers_hist[:,j] + bulge_diskins_hist[:,j])
        ms_tot_bycomp[j,2] = np.sum(bulge_diskins_hist[:,j]) / np.sum(disk_hist[:,j] + bulge_mergers_hist[:,j] + bulge_diskins_hist[:,j])
        ind = np.where(mstars >= 9)
        ms_tot_bycomp_massivegals[j,0] = np.sum(disk_hist[ind,j]) / np.sum(disk_hist[ind,j] + bulge_mergers_hist[ind,j] + bulge_diskins_hist[ind,j])
        ms_tot_bycomp_massivegals[j,1] = np.sum(bulge_mergers_hist[ind,j]) / np.sum(disk_hist[ind,j] + bulge_mergers_hist[ind,j] + bulge_diskins_hist[ind,j])
        ms_tot_bycomp_massivegals[j,2] = np.sum(bulge_diskins_hist[ind,j]) / np.sum(disk_hist[ind,j] + bulge_mergers_hist[ind,j] + bulge_diskins_hist[ind,j])



    return (sfhs_masses, sfr_tot, sfrall, mshs_masses_cum, ms_tot, ms_tot_bycomp, ms_true, ms_tot_bycomp_massivegals, sfhs_comp, mshs_comp)

def plot_cosmic_sfr(plt, outdir, obsdir, sfhs_masses, sfr_tot, sfr_true, mshs_masses, ms_tot, ms_true, LBT, h0):
    
    #load observations
    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    redK11, SFRK11, err_upK11, err_dnK11 = common.load_observation(obsdir, 'Global/SFRD_Karim11.dat', [0,1,2,3])

    hobs = 0.7
    xobs = redK11

    yobs = xobs*0. - 999.
    indx = np.where( SFRK11 > 0)
    yobs[indx] = np.log10(SFRK11[indx] * pow(hobs/h0, 2.0))

    lerr = yobs*0. - 999.
    indx = np.where( (SFRK11-err_dnK11) > 0)
    lerr[indx]  = np.log10(SFRK11[indx] - err_dnK11[indx])

    herr = yobs*0. + 999.
    indx = np.where( (SFRK11+err_upK11) > 0)
    herr[indx]  = np.log10(SFRK11[indx] + err_upK11[indx])

    #Driver (Chabrier IMF), ['Baldry+2012, z<0.06']
    redD17d, redD17u, sfrD17, err1, err2, err3, err4 = common.load_observation(obsdir, 'Global/Driver18_sfr.dat', [0,1,2,3,4,5,6])
    hobs = 0.7
    xobsD17 = (redD17d+redD17u)/2.0
    yobsD17 = sfrD17 + np.log10(hobs/h0)

    errD17 = yobs*0. - 999.
    errD17 = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))

    fig = plt.figure(figsize=(5,4.5))
    ytit = "$\\rm log_{10} (\\rm \\rho_{\\rm SFR}/ M_{\\odot} yr^{-1} cMpc^{-3})$"
    xtit = "lookback time/Gyr"
    xmin, xmax, ymin, ymax = 0, 13, -4, 0
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    #Predicted relation
    cols = ['Khaki','MediumSeaGreen','Turquoise', 'MediumBlue']
    labels = ['$log(M_{\\star})<8.5$', '$8.5<log(M_{\\star})<9.5$','$9.5<log(M_{\\star})<10.5$','$log(M_{\\star})>10.5$']

    for i in range(0,len(cols)):
        ax.plot(LBT, np.log10(sfhs_masses[:,i]), linestyle='solid',color=cols[i], label=labels[i])
    
    print("Will print Cosmic SFR history")
    for a,b,c,d,e,f in zip(LBT, sfr_true[:], sfhs_masses[:,0], sfhs_masses[:,1], sfhs_masses[:,2], sfhs_masses[:,3]):
        print(a,b,c,d,e,f)
 
    ax.plot(LBT, np.log10(sfr_tot[:]), linestyle='solid', color='k')
    ax.plot(LBT, np.log10(sfr_true[:]), linestyle='dashed', color='k')

    ax.errorbar(us.look_back_time(xobs[0:8]), yobs[0:8], yerr=[yobs[0:8]-lerr[0:8],herr[0:8]-yobs[0:8]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D')
    ax.errorbar(us.look_back_time(xobs[9:17]), yobs[9:17], yerr=[yobs[9:17]-lerr[9:17],herr[9:17]-yobs[9:17]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='x')
    ax.errorbar(us.look_back_time(xobsD17), yobsD17, yerr=[errD17,errD17], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

    common.prepare_legend(ax, cols, loc=2)
    common.savefig(outdir, fig, 'cosmic_sfr_forensics.pdf')



    redD17d, redD17u, smdD17, err1, err2, err3, err4 = common.load_observation(obsdir, 'Global/Driver18_smd.dat', [1,2,3,4,5,6,7])
    hobs = 0.7
    xobs = (redD17d+redD17u)/2.0
    yobs = smdD17 + np.log10(hobs/h0)
    err = yobs*0. - 999.
    err = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))

    fig = plt.figure(figsize=(5,4.5))
    ytit = "$\\rm log_{10} (\\rm \\rho_{\\star}/ M_{\\odot} cMpc^{-3})$"
    xtit = "lookback time/Gyr"
    xmin, xmax, ymin, ymax = 0, 13, 4, 9
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    #Predicted relation
    for i in range(0,len(cols)):
        ax.plot(LBT, np.log10(mshs_masses[:,i]), linestyle='solid',color=cols[i], label=labels[i])
 
    ax.plot(LBT, np.log10(ms_tot[:]), linestyle='solid', color='k')
    ax.plot(LBT, np.log10(ms_true[:]), linestyle='dashed', color='k')
    ax.errorbar(us.look_back_time(xobs), yobs, yerr=[err,err], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

    print("Will print Cosmic SMD history")
    for a,b,c,d,e,f in zip(LBT, ms_tot[:], mshs_masses[:,0], mshs_masses[:,1], mshs_masses[:,2], mshs_masses[:,3]):
        print(a,b,c,d,e,f)

    common.prepare_legend(ax, cols, loc=3)
    common.savefig(outdir, fig, 'cosmic_smd_forensics.pdf')

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge')}

    zlist = [0]
    sfh_fields = {'bulges_diskins': ('star_formation_rate_histories'),
                  'bulges_mergers': ('star_formation_rate_histories'),
                  'disks': ('star_formation_rate_histories')}

    fields_global = {'global': ('redshifts', 'mstars', 'sfr_quiescent', 'sfr_burst')}
    # Read data from each subvolume at a time and add it up
    # rather than appending it all together
    for idx, subvol in enumerate(subvols):
        subvol_data = common.read_data(modeldir, redshift_table[0], fields_global, [subvol])
        if idx == 0:
            hdf5_data_hist        = subvol_data
        else:
            for subvol_datum, hdf5_datum in zip(subvol_data[3:], hdf5_data_hist[3:]):
                hdf5_datum += subvol_datum
                #select the most massive black hole from the last list item

    # Also make sure that the total volume takes into account the number of subvolumes read
    hdf5_data_hist[1] = hdf5_data_hist[1] * len(subvols)
    h0, redshifts = hdf5_data_hist[0], hdf5_data_hist[2]

    hdf5_data = common.read_data(modeldir, 199, fields, subvols)
    sfh, delta_t, LBT = common.read_sfh(modeldir, 199, sfh_fields, subvols)
    (sfhs_masses, sfr_tot, sfr_true, mshs_masses, ms_tot, ms_tot_bycomp, ms_true, ms_tot_bycomp_massivegals, sfhs_comp, mshs_comp) = prepare_data(hdf5_data, sfh, hdf5_data_hist, LBT, delta_t)

    print("#SFRD Forensics")
    print("#LBT/Gyr SFRD_disk[Msun/yr/Mpc^3] SFRD_bulge_mergers[Msun/yr/Mpc^3] SFRD_bulge_diskins[Msun/yr/Mpc^3]")
    for a,b,c,d in zip(LBT, sfhs_comp[:,0], sfhs_comp[:,1], sfhs_comp[:,2]):
        print(a,b,c,d)

    print("#SMD Forensics")
    print("#LBT/Gyr SMD_disk[Msun/Mpc^3] SMD_bulge_mergers[Msun/Mpc^3] SMD_bulge_diskins[Msun/Mpc^3]")
    for a,b,c,d in zip(LBT, mshs_comp[:,0], mshs_comp[:,1], mshs_comp[:,2]):
        print(a,b,c,d)
       

    plot_cosmic_sfr(plt, outdir, obsdir, sfhs_masses, sfr_tot, sfr_true, mshs_masses, ms_tot, ms_true, LBT, h0)


if __name__ == '__main__':
    main(*common.parse_args())
