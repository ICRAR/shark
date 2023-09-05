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
import os

import common


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

def plot_uv_lf_evo(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust, nbands):

    volcorr = 3.0*np.log10(h0)
    xlf_obs  = xlf
 
    xtit="$\\rm 1500\\AA\, mag\, (AB)$"
    ytit="$\\rm log_{10}(\Phi/{\\rm dex^{-1}} {\\rm Mpc}^{-3})$"

    xmin, xmax, ymin, ymax = -25, -15, -6, -1
    xleg = xmin + 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    fig = plt.figure(figsize=(5,14))

    subplots = (511, 512, 513, 514, 515)
    idx = (0, 1, 2, 3, 4)
    zs  = (0, 1, 2, 3, 4)
    band = 0 #28
    labels= ('z=3', 'z=4', 'z=6', 'z=8', 'z=10')
  
    corrm_obs = -5.0*np.log10(h0/0.7) 
    corry_obs = 3.0*np.log10(h0/0.7)
    for subplot, idx, z in zip(subplots, idx, zs):

        ax = fig.add_subplot(subplot)
        ytitplot = ytit
        if (idx == 4):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        if(idx == 0):
           file = obsdir+'/lf/lf1700_z3_sawicki06.data'
           lm,p,dp = np.loadtxt(file,usecols=[0,2,3],unpack=True)
           indx = np.where(p > 0)
           yobs = np.log10(p[indx])
           ydn  = np.log10(p[indx]-dp[indx])
           yup  = np.log10(p[indx]+dp[indx])
           ax.errorbar(lm[indx]+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Sawicki+2006")

           file = obsdir+'/lf/lf1700_z3_reddy09.data'
           lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
           indx = np.where(p > 0)
           yobs = np.log10(p[indx])
           ydn  = np.log10(p[indx]-dp[indx])
           yup  = np.log10(p[indx]+dp[indx])
           ax.errorbar(lm[indx]+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Reddy+2009")

        if(idx == 1):
           file = obsdir+'/lf/lf1500_z4_adams19.data'
           lmA19z4,Ap4,Adp4 = np.loadtxt(file,usecols=[0, 1, 2],unpack=True)
           yobs = np.log10(Ap4*1e-4)
           ydn  = np.log10(Ap4*1e-4 - Adp4*1e-4)
           yup  = np.log10(Ap4*1e-4 + Adp4*1e-4)
           ax.errorbar(lmA19z4+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='d',label="Adams+2019")

           file = obsdir+'/lf/lf1600_z4-10_Bouwens2015.data'
           lmB15z4,p4,dp4,lmB15z6,p6,dp6,lmB15z8,p8,dp8 = np.loadtxt(file,usecols=[0, 1, 2, 6, 7, 8, 12, 13, 14],unpack=True)
           yobs = np.log10(p4)
           ydn  = np.log10(p4-dp4)
           yup  = np.log10(p4+dp4)
           ax.errorbar(lmB15z4+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^',label="Bouwens+2015")

           file = obsdir+'/lf/lf1500_z4-8_Finkelstein2015.data'
           lmF15,pF4,dpuF4,dpdF4,pF6,dpuF6,dpdF6,pF8,dpuF8,dpdF8 = np.loadtxt(file,usecols=[0,1, 2, 3, 7, 8, 9, 13, 14, 15],unpack=True)
           yobs = np.log10(pF4*1e-3)
           ydn  = np.log10(pF4*1e-3-dpdF4*1e-3)
           yup  = np.log10(pF4*1e-3+dpuF4*1e-3)
           ax.errorbar(lmF15+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='v',label="Finkelstein+2015")

        if(idx == 2):
           yobs = np.log10(p6)
           ydn  = np.log10(p6-dp6)
           yup  = np.log10(p6+dp6)
           ax.errorbar(lmB15z6+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^')

           yobs = np.log10(pF6*1e-3)
           ydn  = np.log10(pF6*1e-3-dpdF6*1e-3)
           yup  = np.log10(pF6*1e-3+dpuF6*1e-3)
           ax.errorbar(lmF15+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='v')

        if(idx == 3):
           yobs = np.log10(p8)
           ydn  = np.log10(p8-dp8)
           yup  = np.log10(p8+dp8)
           ax.errorbar(lmB15z8+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^')

           yobs = np.log10(pF8*1e-3)
           ydn  = np.log10(pF8*1e-3-dpdF8*1e-3)
           yup  = np.log10(pF8*1e-3+dpuF8*1e-3)
           ax.errorbar(lmF15+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='v')

           file = obsdir+'/lf/lf1500_z8_adams23.data'
           lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
           ax.errorbar(lm+corrm_obs, np.log10(p * 1e-5)+corry_obs, yerr=[np.log10(p)-np.log10(dp),np.log10(p + dp) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D', label='Adams+2023')

        if(idx == 4):
           file = obsdir+'/lf/lf1500_z10_oesch2018.data'
           lm,p,dpu,dpd = np.loadtxt(file,usecols=[0,1, 2, 3],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[p-dpu,dpd-p], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*', label='Oesch+2018')

           file = obsdir+'/lf/lf1500_z10_adams23.data'
           lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
           ax.errorbar(lm+corrm_obs, np.log10(p * 1e-5)+corry_obs, yerr=[np.log10(p)-np.log10(dp),np.log10(p + dp) - np.log10(p)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D')

        #Predicted LF
        ind = np.where(LFs_dust[z,4,band,:] < 0.)
        y = LFs_dust[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
        ind = np.where(LFs_nodust[z,4,band,:] < 0.)
        y = LFs_nodust[z,4,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)
        if(idx == 1):
           for a,b,c in zip(xlf_obs,LFs_dust[z,4,band,:],LFs_nodust[z,4,band,:]):
               print (a, b+volcorr-np.log10(dm),c+volcorr-np.log10(dm))

        ind = np.where(LFs_dust[z,3,band,:] < 0.)
        y = LFs_dust[z,3,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
        ind = np.where(LFs_dust[z,2,band,:] < 0.)
        y = LFs_dust[z,2,band,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
        if ((idx == 0) or (idx == 1) or (idx ==4) or (idx == 3)):
            common.prepare_legend(ax, ['grey','grey','grey'], loc=4)

    plt.tight_layout()
    common.savefig(outdir, fig, "UV_luminosity_function_evolution.pdf")

def prepare_data(hdf5_data, phot_data, phot_data_nod, LFs_dust, LFs_nodust, index, nbands):
   
    #star_formation_histories and SharkSED have the same number of galaxies in the same order, and so we can safely assume that to be the case.
    #to select the same galaxies in galaxies.hdf5 we need to ask for all of those that have a stellar mass > 0, and then assume that they are in the same order.

    (h0, _, mdisk, mbulge, mhalo, mshalo, typeg, age, 
     sfr_disk, sfr_burst, id_gal) = hdf5_data
   
    #(bulge_diskins_hist, bulge_mergers_hist, disk_hist) = sfh

    #components:
    #(len(my_data), 2, 2, 5, nbands)
    #0: disk instability bulge
    #1: galaxy merger bulge
    #2: total bulge
    #3: disk
    #4: total

    ind = np.where(mdisk + mbulge > 0)
    SEDs_dust = np.zeros(shape = (len(mdisk[ind]), 5, nbands))
    SEDs_nodust = np.zeros(shape = (len(mdisk[ind]), 5, nbands))

    p = 0
    for c in range(0,5):
        indust = phot_data[p]
        innodust = phot_data_nod[p]
        for i in range(0,nbands):
            SEDs_dust[:,c,i] = indust[i,:]
            SEDs_nodust[:,c,i] = innodust[i,:]
        p = p + 1

    for i in range(0,nbands):
        for c in range(0,5):
            #calculate LF with bands with dust
            ind = np.where(SEDs_dust[:,c,i] < -1)
            H, bins_edges = np.histogram(SEDs_dust[ind,c,i],bins=np.append(mbins,mupp))
            LFs_dust[index,c,i,:] = LFs_dust[index,c,i,:] + H

            #calculate LF of intrinsic bands 
            ind = np.where(SEDs_nodust[:,c,i] < -1)
            H, bins_edges = np.histogram(SEDs_nodust[ind,c,i],bins=np.append(mbins,mupp))
            LFs_nodust[index,c,i,:] = LFs_nodust[index,c,i,:] + H


def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age', 
                           'sfr_disk', 'sfr_burst', 'id_galaxy')}

    #sfh_fields = {'bulges_diskins': ('star_formation_rate_histories'),
    #              'bulges_mergers': ('star_formation_rate_histories'),
    #              'disks': ('star_formation_rate_histories')}

    Variable_Ext = True

    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}
    fields_sed_nod = {'SED/ab_nodust': ('bulge_d','bulge_m','bulge_t','disk','total')}

    z = (3.0, 4.0, 6.0, 8.0, 10.0)  #(8.0, 9, 10.5, 12.5) #, 1.0, 1.5, 2.0)
    snapshots = redshift_table[z]

    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5" 
    # Create histogram
    for index, snapshot in enumerate(snapshots):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        if(Variable_Ext == False):
           seds = common.read_photometry_data(model_dir, snapshot, fields_sed, subvols)
           seds_nod = common.read_photometry_data(model_dir, snapshot, fields_sed_nod, subvols)
        else:
           seds = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed)
           seds_nod = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_nod, subvols, file_hdf5_sed)

        nbands = len(seds[0]) 

        print(nbands)
        if(index == 0):
            LFs_dust     = np.zeros(shape = (len(z), 5, nbands, len(mbins)))
            LFs_nodust   = np.zeros(shape = (len(z), 5, nbands, len(mbins)))

        prepare_data(hdf5_data, seds, seds_nod, LFs_dust, LFs_nodust, index, nbands)

        h0, volh = hdf5_data[0], hdf5_data[1]
        if(volh > 0.):
            LFs_dust[index,:]   = LFs_dust[index,:]/volh
            LFs_nodust[index,:] = LFs_nodust[index,:]/volh

    print ("number of bands %d" % (nbands,))
    # Take logs
    ind = np.where(LFs_dust > 0.)
    LFs_dust[ind] = np.log10(LFs_dust[ind])

    ind = np.where(LFs_nodust > 0.)
    LFs_nodust[ind] = np.log10(LFs_nodust[ind])

    if(Variable_Ext):
       outdir = os.path.join(outdir, 'eagle-rr14')

    volcorr = 3.0*np.log10(h0)

    band = 0
    for a,b,c,d,e in zip(xlf, LFs_dust[0,4,band,:], LFs_dust[1,4,band,:],LFs_dust[2,4,band,:],LFs_dust[3,4,band,:]):
        print(a,b+volcorr-np.log10(dm),c+volcorr-np.log10(dm),d+volcorr-np.log10(dm),e+volcorr-np.log10(dm))
    plot_uv_lf_evo(plt, outdir, obsdir, h0, LFs_dust, LFs_nodust, nbands)

if __name__ == '__main__':
    main(*common.parse_args())
