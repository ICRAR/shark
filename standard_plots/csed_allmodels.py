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

def plot_csed(plt, outdir, obsdir, h0, CSED, CSED_nodust, CSED2, CSED3, CSED4, nbands):

    #wavelength in angstroms.
    file = obsdir+'/Models/Shark_SED_bands.dat'
    lambda_bands = np.loadtxt(file,usecols=[0],unpack=True)
    freq_bands   = c_light / (lambda_bands * 1e-10) #in Hz
    lambda_bands = np.log10(lambda_bands)

    xtit="$\\rm log_{10}(\lambda/Ang\, (rest-frame))$"
    ytit="$\\rm log_{10}(\\nu \\epsilon_{\\rm int}/ h\,W\, Mpc^{-3})$"

    fig = plt.figure(figsize=(6,13))

    subplots = (411, 412, 413, 414)
    idx = (0, 1, 2, 3)
    labels= ('z=0.25', 'z=1', 'z=3', 'z=8')
    obs = ('z0p25', 'z0p95')
    colors = ('Indigo','purple','Navy','MediumBlue','Green','MediumAquamarine','LightGreen','YellowGreen','Gold','Orange','Coral','OrangeRed','red','DarkRed','FireBrick','Crimson','IndianRed','LightCoral','Maroon','brown','Sienna','SaddleBrown','Chocolate','Peru','DarkGoldenrod','Goldenrod','SandyBrown')

    for subplot, idx in zip(subplots, idx):
        xmin, xmax, ymin, ymax = 3.0, 7.0, 32, 36
        xleg = xmin + 0.1 * (xmax-xmin)
        yleg = ymin + 0.1 * (ymax-ymin)

        ax = fig.add_subplot(subplot)
        if (idx == 3):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytit, locators=(1, 1, 1, 1))
        ax.text(xleg,yleg, labels[idx], fontsize=12)

        if(idx <= 1):
           #plot observations
           file = obsdir+'/lf/CSED/CSED_Andrews17_'+obs[idx]+'.dat'
           lw,p1,p2 = np.loadtxt(file,usecols=[0,1,2],unpack=True)
           lw = np.log10(lw * 1e10) #in angstroms 
           p1 = np.log10(p1)
           p2 = np.log10(p2)
      
           xobs = np.zeros(shape = 2)
           yobs = np.zeros(shape = 2)
           for xi,ymin,ymax in zip(lw,p1,p2):
               xobs[0] = xi
               xobs[1] = xi
               yobs[0] = ymin
               yobs[1] = ymax
               ax.plot(xobs,yobs, color='grey', linestyle='solid',linewidth=5)  

        #plot model
        ax.plot(lambda_bands,np.log10(CSED[idx,4,:]*freq_bands)-np.log10(h0), 'Indigo', linewidth=3, alpha=0.5, label='$\\rm EAGLE-\\tau\, RR14-steep$')
        ax.plot(lambda_bands,np.log10(CSED_nodust[idx,4,:]*freq_bands)-np.log10(h0), 'Indigo', linewidth=1, label='intrinsic')
        ax.plot(lambda_bands,np.log10(CSED2[idx,4,:]*freq_bands)-np.log10(h0), 'Indigo', linewidth=3, linestyle='dotted', label='$\\rm EAGLE-\\tau\, RR14$')
        ax.plot(lambda_bands,np.log10(CSED3[idx,4,:]*freq_bands)-np.log10(h0), 'Indigo', linewidth=3, alpha=0.4, linestyle='dashed', label='$\\rm EAGLE-\\tau\,f_{\\rm dust}\, const$')
        ax.plot(lambda_bands,np.log10(CSED4[idx,4,:]*freq_bands)-np.log10(h0), 'Indigo', linewidth=3, alpha=0.3, linestyle='dashdot', label='CF00')

        if (idx == 0):
            common.prepare_legend(ax, ['Indigo','Indigo','Indigo','Indigo','Indigo'], bbox_to_anchor=[0.2, 1], handlelength=5)

    common.savefig(outdir, fig, "CSED_Shark_allmodels.pdf")

def prepare_data(hdf5_data, phot_data, phot_data_nodust, seds2, seds3, seds4, CSED, CSED_nodust, CSED2, CSED3, CSED4, nbands, index):
   
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

    z = (0.25, 1, 3.0, 8.0)
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

        prepare_data(hdf5_data, seds, seds_nodust, seds2, seds3, seds4, CSED, CSED_nodust, CSED2, CSED3, CSED4, nbands, index)

        h0, volh = hdf5_data[0], hdf5_data[1]
        if(volh > 0.):
            CSED[index,:]   = CSED[index,:] / volh * pow(h0,3.0)
            CSED_nodust[index,:]   = CSED_nodust[index,:] / volh * pow(h0,3.0)
            CSED2[index,:]   = CSED2[index,:] / volh * pow(h0,3.0)
            CSED3[index,:]   = CSED3[index,:] / volh * pow(h0,3.0)
            CSED4[index,:]   = CSED4[index,:] / volh * pow(h0,3.0)

    if(Variable_Ext):
       outdir = os.path.join(outdir, 'eagle-rr14-steep')

    # Take logs
    plot_csed(plt, outdir, obsdir, h0, CSED, CSED_nodust, CSED2, CSED3, CSED4, nbands)
 
if __name__ == '__main__':
    main(*common.parse_args())
