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

# colour distribution initialization
clow  = -0.1
cupp  = 3.5
dc    = 0.075
cbins = np.arange(clow,cupp,dc)
xc    = cbins + dc/2.0

magbins = [-17.13,-17.88,-18.63,-19.38,-20.13,-20.88,-21.63]


def plot_uvslope(plt, outdir, obsdir, h0, z, beta_uv, beta_uv2, beta_uv3, beta_uv4,
                beta_uv_mag, beta_uv2_mag,  beta_uv3_mag, beta_uv4_mag):

    xtit="$\\rm redshift$"
    ytit="$\\beta_{\\rm UV,CSED}$"

    fig = plt.figure(figsize=(5,4))
    plt.subplots_adjust(left=0.15, bottom=0.17)

    xmin, xmax, ymin, ymax = 0,10,-2.6,0.1

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

    ax.plot(z, beta_uv, 'Indigo', linewidth=3, alpha=0.5, label='$\\rm EAGLE-\\tau\, RR14-steep$')
    ax.plot(z, beta_uv2, 'Indigo', linewidth=3, linestyle='dotted', label='$\\rm EAGLE-\\tau\, RR14$')
    ax.plot(z, beta_uv3, 'Indigo', linewidth=3, alpha=0.4, linestyle='dashed', label='$\\rm EAGLE-\\tau\,f_{\\rm dust}\, const$')
    ax.plot(z, beta_uv4, 'Indigo', linewidth=3, alpha=0.3, linestyle='dashdot', label='CF00')

    x = [3,3]
    ye = [0.25,0.25]
    y = [-1.55,-1.55]
    ax.errorbar(x, y, yerr=[ye,ye], ls='None', mfc='None', ecolor = 'k', mec='k',marker='o', label='Davies+2013')

    common.prepare_legend(ax, ['Indigo','Indigo','Indigo','Indigo','k'], loc='upper right', handlelength=2.5)
    common.savefig(outdir, fig, "UV_Slope_CSED_Shark_allmodels.pdf")

    #UV slope at -19.5
    xtit="$\\rm redshift$"
    ytit="$\\beta_{\\rm UV,-19.5}$"

    fig = plt.figure(figsize=(5,4))
    plt.subplots_adjust(left=0.15, bottom=0.17)

    xmin, xmax, ymin, ymax = 0,10,-3,0.1

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

    ax.plot(z, beta_uv_mag, 'Indigo', linewidth=3, alpha=0.5)
    ax.plot(z, beta_uv2_mag, 'Indigo', linewidth=3, linestyle='dotted')
    ax.plot(z, beta_uv3_mag, 'Indigo', linewidth=3, alpha=0.4, linestyle='dashed')
    ax.plot(z, beta_uv4_mag, 'Indigo', linewidth=3, alpha=0.3, linestyle='dashdot')

    file = obsdir+'/lf/beta_uv_19p5_Bouwens14.data'
    lm,p,dpu,dpd = np.loadtxt(file,usecols=[0,1, 2, 3],unpack=True)
    ax.errorbar(lm, p, yerr=[dpu+dpd,dpu+dpd], ls='None', mfc='None', ecolor = 'k', mec='k',marker='*', label='Bouwens+2014')

    common.prepare_legend(ax, ['k'], loc='upper right')
    common.savefig(outdir, fig, "UV_Slope_m19p5_Shark_allmodels.pdf")

def prepare_data(hdf5_data, phot_data, phot_data_nodust, seds2, seds3, seds4, CSED, CSED_nodust, CSED2, CSED3, CSED4, 
                 nbands, index,  beta_uv, beta_uv2, beta_uv3, beta_uv4, 
                 beta_uv_mag, beta_uv2_mag,  beta_uv3_mag, beta_uv4_mag, obsdir):
   
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
    SEDs_dust = np.zeros(shape = (len(mdisk[ind]), 5, nbands-1))
    SEDs_nodust = np.zeros(shape = (len(mdisk[ind]), 5, nbands-1))
    SEDs_dust2 = np.zeros(shape = (len(mdisk[ind]), 5, nbands-1))
    SEDs_dust3 = np.zeros(shape = (len(mdisk[ind]), 5, nbands-1))
    SEDs_dust4 = np.zeros(shape = (len(mdisk[ind]), 5, nbands-1))

    p = 0
    for c in range(0,5):
        indust = phot_data[p]
        nodust = phot_data_nodust[p]
        indust2 = seds2[p]
        indust3 = seds3[p]
        indust4 = seds4[p]

        for i in range(0,nbands-1):
            SEDs_dust[:,c,i] = indust[i,:]
            SEDs_nodust[:,c,i] = nodust[i,:]
            SEDs_dust2[:,c,i] = indust2[i,:]
            SEDs_dust3[:,c,i] = indust3[i,:]
            SEDs_dust4[:,c,i] = indust4[i,:]

        p = p + 1

    totuv = np.zeros(shape = (3,4))

    for i in range(0,nbands-1):
        for c in range(0,5):
            #calculate LF with bands with dust
            ind = np.where((SEDs_dust[:,c,i] < -1) & (SEDs_dust[:,c,i] > -50))
            #W Hz-1
            total = np.sum(pow(10.0,(SEDs_dust[ind,c,i]+48.6)/(-2.5))) * (4.0 * PI * pow(10.0*3.086e18,2.0)) * 1e-7 
            CSED[index,c,i] = total
            #no dust 
            ind = np.where((SEDs_nodust[:,c,i] < -1) & (SEDs_nodust[:,c,i] > -50))
            total = np.sum(pow(10.0,(SEDs_nodust[ind,c,i]+48.6)/(-2.5))) * (4.0 * PI * pow(10.0*3.086e18,2.0)) * 1e-7 
            CSED_nodust[index,c,i] = total 
            #different extinction models
            ind = np.where((SEDs_dust2[:,c,i] < -1) & (SEDs_dust2[:,c,i] > -50))
            total = np.sum(pow(10.0,(SEDs_dust2[ind,c,i]+48.6)/(-2.5))) * (4.0 * PI * pow(10.0*3.086e18,2.0)) * 1e-7 
            CSED2[index,c,i] = total
            ind = np.where((SEDs_dust3[:,c,i] < -1) & (SEDs_dust3[:,c,i] > -50))
            total = np.sum(pow(10.0,(SEDs_dust3[ind,c,i]+48.6)/(-2.5))) * (4.0 * PI * pow(10.0*3.086e18,2.0)) * 1e-7 
            CSED3[index,c,i] = total
            ind = np.where((SEDs_dust4[:,c,i] < -1) & (SEDs_dust4[:,c,i] > -50))
            total = np.sum(pow(10.0,(SEDs_dust4[ind,c,i]+48.6)/(-2.5))) * (4.0 * PI * pow(10.0*3.086e18,2.0)) * 1e-7 
            CSED4[index,c,i] = total

            #compute required quantities for UV slope of galaxies of mag ~-19.5
            if((c == 4) & (i <= 2)):
               ind = np.where((SEDs_dust[:,c,i] < -19.3) & (SEDs_dust[:,c,i] > -19.7))  
               totuv[i,0] = np.sum(pow(10.0,(SEDs_dust[ind,c,i]+48.6)/(-2.5))) * (4.0 * PI * pow(10.0*3.086e18,2.0)) * 1e-7
               ind = np.where((SEDs_dust2[:,c,i] < -19.3) & (SEDs_dust2[:,c,i] > -19.7))  
               totuv[i,1] = np.sum(pow(10.0,(SEDs_dust2[ind,c,i]+48.6)/(-2.5))) * (4.0 * PI * pow(10.0*3.086e18,2.0)) * 1e-7
               ind = np.where((SEDs_dust3[:,c,i] < -19.3) & (SEDs_dust3[:,c,i] > -19.7))  
               totuv[i,2] = np.sum(pow(10.0,(SEDs_dust3[ind,c,i]+48.6)/(-2.5))) * (4.0 * PI * pow(10.0*3.086e18,2.0)) * 1e-7
               ind = np.where((SEDs_dust4[:,c,i] < -19.3) & (SEDs_dust4[:,c,i] > -19.7))  
               totuv[i,3] = np.sum(pow(10.0,(SEDs_dust4[ind,c,i]+48.6)/(-2.5))) * (4.0 * PI * pow(10.0*3.086e18,2.0)) * 1e-7
             

    file = obsdir+'/Models/Shark_SED_bands.dat'
    lambda_bands = np.loadtxt(file,usecols=[0],unpack=True)
    freq_bands   = c_light / (lambda_bands * 1e-10) #in Hz
    lambda_bands = np.log10(lambda_bands)

    ind = np.where((lambda_bands > 3)  & (lambda_bands < 3.56))
    y = CSED[index,4,ind]
    (uv_fit_slope, uv_fit_offs) = np.polyfit(lambda_bands[ind],np.log10(y[0]),1)
    beta_uv[index] = uv_fit_slope-2
    y = CSED2[index,4,ind]
    (uv_fit_slope, uv_fit_offs) = np.polyfit(lambda_bands[ind],np.log10(y[0]),1)
    beta_uv2[index] = uv_fit_slope-2
    y = CSED3[index,4,ind]
    (uv_fit_slope, uv_fit_offs) = np.polyfit(lambda_bands[ind],np.log10(y[0]),1)
    beta_uv3[index] = uv_fit_slope-2
    y = CSED4[index,4,ind]
    (uv_fit_slope, uv_fit_offs) = np.polyfit(lambda_bands[ind],np.log10(y[0]),1)
    beta_uv4[index] = uv_fit_slope-2

    #UV slope at -19.5
    (uv_fit_slope, uv_fit_offs) = np.polyfit(lambda_bands[ind],np.log10(totuv[:,0]),1)
    beta_uv_mag[index] = uv_fit_slope-2
    (uv_fit_slope, uv_fit_offs) = np.polyfit(lambda_bands[ind],np.log10(totuv[:,1]),1)
    beta_uv2_mag[index] = uv_fit_slope-2
    (uv_fit_slope, uv_fit_offs) = np.polyfit(lambda_bands[ind],np.log10(totuv[:,2]),1)
    beta_uv3_mag[index] = uv_fit_slope-2
    (uv_fit_slope, uv_fit_offs) = np.polyfit(lambda_bands[ind],np.log10(totuv[:,3]),1)
    beta_uv4_mag[index] = uv_fit_slope-2

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    Variable_Ext = True 
    file_hdf5_sed = "Shark-SED-eagle-rr14-steep.hdf5"
    file_hdf5_sed2 = "Shark-SED-eagle-rr14.hdf5"
    file_hdf5_sed3 = "Shark-SED-eagle-const.hdf5" 
    file_hdf5_sed4 = "Shark-SED.hdf5"

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age', 
                           'sfr_disk', 'sfr_burst', 'id_galaxy')}
    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}
    fields_sed_nodust = {'SED/ab_nodust': ('bulge_d','bulge_m','bulge_t','disk','total'),}

    z = (0, 0.25, 0.5, 1, 1, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0)
    snapshots = redshift_table[z]

    # Create histogram
    for index, snapshot in enumerate(snapshots):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        #sfh, delta_t, LBT = common.read_sfh(model_dir, snapshot, sfh_fields, subvols)
        if(Variable_Ext == False):
           seds = common.read_photometry_data(model_dir, snapshot, fields_sed, subvols)
           seds_nodust = common.read_photometry_data(model_dir, snapshot, fields_sed_nodust, subvols)
        else:
           seds = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed)
           seds_nodust = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_nodust, subvols, file_hdf5_sed)
           seds2 = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed2)
           seds3 = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed3)
           seds4 = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed4)

        nbands = len(seds[0]) 

        if(index == 0):
            CSED = np.zeros(shape = (len(z), 5, nbands-1))
            CSED_nodust = np.zeros(shape = (len(z), 5, nbands-1))
            CSED2 = np.zeros(shape = (len(z), 5, nbands-1))
            CSED3 = np.zeros(shape = (len(z), 5, nbands-1))
            CSED4 = np.zeros(shape = (len(z), 5, nbands-1))
            beta_uv = np.zeros(shape = (len(z)))
            beta_uv2 = np.zeros(shape = (len(z)))
            beta_uv3 = np.zeros(shape = (len(z)))
            beta_uv4 = np.zeros(shape = (len(z)))
            beta_uv_mag = np.zeros(shape = (len(z)))
            beta_uv2_mag = np.zeros(shape = (len(z)))
            beta_uv3_mag = np.zeros(shape = (len(z)))
            beta_uv4_mag = np.zeros(shape = (len(z)))

        prepare_data(hdf5_data, seds, seds_nodust, seds2, seds3, seds4, CSED, CSED_nodust, CSED2, CSED3, CSED4, 
                     nbands, index, beta_uv, beta_uv2, beta_uv3, beta_uv4, beta_uv_mag, 
                     beta_uv2_mag,  beta_uv3_mag, beta_uv4_mag, obsdir)

        h0, volh = hdf5_data[0], hdf5_data[1]
        if(volh > 0.):
            CSED[index,:]   = CSED[index,:] / volh * pow(h0,3.0)
            CSED_nodust[index,:]   = CSED_nodust[index,:] / volh * pow(h0,3.0)
            CSED2[index,:]   = CSED2[index,:] / volh * pow(h0,3.0)
            CSED3[index,:]   = CSED3[index,:] / volh * pow(h0,3.0)
            CSED4[index,:]   = CSED4[index,:] / volh * pow(h0,3.0)

    if(Variable_Ext):
       outdir = os.path.join(outdir, 'eagle-rr14-steep')

    plot_uvslope(plt, outdir, obsdir, h0, z, beta_uv, beta_uv2, beta_uv3, beta_uv4,
                beta_uv_mag, beta_uv2_mag,  beta_uv3_mag, beta_uv4_mag)
 
if __name__ == '__main__':
    main(*common.parse_args())
