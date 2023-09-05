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

def prepare_data(hdf5_data, phot_data, index, nbands, redshift):
    #star_formation_histories and SharkSED have the same number of galaxies in the same order, and so we can safely assume that to be the case.
    #to select the same galaxies in galaxies.hdf5 we need to ask for all of those that have a stellar mass > 0, and then assume that they are in the same order.

    (h0, _, mdisk, mbulge, mhalo, mshalo, typeg, age, 
     sfr_disk, sfr_burst, id_gal) = hdf5_data

    seds_tot =  phot_data[4] #total absolute magnitudes with dust
   
    #(bulge_diskins_hist, bulge_mergers_hist, disk_hist) = sfh

    #components:
    #(len(my_data), 2, 2, 5, nbands)
    #0: disk instability bulge
    #1: galaxy merger bulge
    #2: total bulge
    #3: disk
    #4: total
    band1 = 1
    band2 = 2 
    band4 = 8

    Vdust = seds_tot[4,:] + 0.4424 * (seds_tot[3,:] - seds_tot[4,:]) + 0.028

    ind = np.where(mdisk + mbulge > 0)
    mstar = np.log10((mdisk[ind] + mbulge[ind])/h0)
    sfr = (sfr_disk[ind] + sfr_burst[ind])/h0/1e9
    ind = np.where(sfr <= 1e-4)    
    sfr[ind] = 1e-4

    ind = np.where(mstar > 9)
    mstarin = mstar[ind]
    sfrin = sfr[ind]
    m1 = seds_tot[band1, ind]
    m2 = seds_tot[band2, ind]
    m3 = Vdust[ind]
    m4 = seds_tot[band4, ind]
    m1 = m1[0]
    m2 = m2[0]
    m4 = m4[0]

    print(len(mstarin))
    ind = np.where((mstar > 9) & (np.log10(sfr)-mstar < -13))
    print(len(mstar[ind]))
 
    print(max(mstarin))
    print("#Mstar[logMsun] > 9 at redshift ", redshift)
    print("#Mstar[logMsun] SFR[Msun/yr] NUV_abs u_abs V_abs J_abs")
    for a,b,c,d,e,f in zip(mstarin, sfrin, m1, m2, m3, m4):
        print(a,b,c,d,e,f)
 

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age', 
                           'sfr_disk', 'sfr_burst', 'id_galaxy')}

    #Bands information:
    #(0): "FUV_GALEX", "NUV_GALEX", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS",
    #(6): "z_SDSS", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1_WISE",
    #(12): "I1_Spitzer", "I2_Spitzer", "W2_WISE", "I3_Spitzer", "I4_Spitzer",
    #(17): "W3_WISE", "W4_WISE", "P70_Herschel", "P100_Herschel",
    #(21): "P160_Herschel", "S250_Herschel", "S350_Herschel", "S450_JCMT",
    #(25): "S500_Herschel", "S850_JCMT", "Band9_ALMA", "Band8_ALMA",
    #(29): "Band7_ALMA", "Band6_ALMA", "Band5_ALMA", "Band4_ALMA"

    #sfh_fields = {'bulges_diskins': ('star_formation_rate_histories'),
    #              'bulges_mergers': ('star_formation_rate_histories'),
    #              'disks': ('star_formation_rate_histories')}

    Variable_Ext = True

    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}

    z = [0] #, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0) #, 1.0, 1.5, 2.0)
    snapshots = redshift_table[z]

    file_hdf5_sed = "Shark-SED-eagle-rr14-steep.hdf5"
    file_hdf5_sed2 = "Shark-SED-eagle-rr14.hdf5"
    file_hdf5_sed3 = "Shark-SED-eagle-const.hdf5"
    file_hdf5_sed4 = "Shark-SED.hdf5"

    # Create histogram
    for index, snapshot in enumerate(snapshots):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        if(Variable_Ext == False):
           seds = common.read_photometry_data(model_dir, snapshot, fields_sed, subvols)
        else:
           seds = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed2)

        nbands = len(seds[0]) 

        prepare_data(hdf5_data, seds, index, nbands, z[index])


if __name__ == '__main__':
    main(*common.parse_args())
