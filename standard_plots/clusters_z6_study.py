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

import functools

import numpy as np
import h5py

import common
import utilities_statistics as us

##################################
# Constants
mlow = 7.0
mupp = 13.0
dm = 0.5
mbins = np.arange(mlow, mupp, dm)
xmf = mbins + dm/2.0

#rbins = np.array([1, 2.730045591676135, 5]) #Mpc/h
rbins = np.array([1, 2.8315841879187973, 5]) #Mpc/h
zdepth = 40.30959350543804 #Mpc/h

GyrtoYr  = 1e9
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun
PI       = 3.1416

def add_observations_to_plot(obsdir, fname, ax, marker, label, color='k', err_absolute=False):
    fname = '%s/Gas/%s' % (obsdir, fname)
    x, y, yerr_down, yerr_up = common.load_observation(obsdir, fname, (0, 1, 2, 3))
    common.errorbars(ax, x, y, yerr_down, yerr_up, color, marker, label=label, err_absolute=err_absolute)

def prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit):
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    xleg = xmax - 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)
    #ax.text(xleg, yleg, 'z=0')

def prepare_data(hdf5_data, seds_ap, seds_ab, index, zlist):


    # Unpack data
    (h0, _, typeg, mdisk, mbulge, _, _, mHI, mH2, mgas,
     mHI_bulge, mH2_bulge, mgas_bulge, mvir, sfrd, sfrb, 
     x, y, z, vvir, vx, vy, vz, id_halo_tree) = hdf5_data
    ap_mags = seds_ap[0]
    ab_mags = seds_ab[0]

    ind = np.where( (mdisk +  mbulge) > 0)
    mvir = mvir[ind]
    vvir = vvir[ind]
    mH2 = mH2[ind]
    mH2_bulge = mH2_bulge[ind]
    x = x[ind]
    y = y[ind]
    z = z[ind]
    vx = vx[ind]
    vy = vy[ind]
    vz = vz[ind]
    mdisk = mdisk[ind]
    mbulge = mbulge[ind]
    sfrd = sfrd[ind]
    sfrb = sfrb[ind]
    typeg = typeg[ind]
    id_halo_tree = id_halo_tree[ind]

    XH = 0.72
    h0log = np.log10(float(h0))

    rvir  = G * mvir / pow(vvir,2.0) / h0

    mstar_tot = (mdisk + mbulge) / h0
    sfr_tot = (sfrd + sfrb) / h0 / GyrtoYr

    #find most massive halos
    indcen = np.where((mvir/h0 > 1e11) & (typeg == 0))
    mvir_in = mvir[indcen]/h0
    sortedmass = np.argsort(1.0/mvir_in)
    mass_ordered = mvir_in[sortedmass]
    mass_tresh  = mass_ordered[9]

    print(max(mvir/h0), mass_ordered)
    indcen = np.where((mvir/h0 >= mass_tresh) & (typeg == 0))
    print(mvir[indcen]/h0)
    x_cen = x[indcen]
    y_cen = y[indcen]
    z_cen = z[indcen]
    vx_cen = vx[indcen]
    vy_cen = vy[indcen]
    vz_cen = vz[indcen]
    mvir_cen = mvir[indcen]/h0
    rvir_cen = rvir[indcen]

    #Bands information:
    #(0): "FUV_GALEX", "NUV_GALEX", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS",
    #(6): "z_SDSS", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1_WISE",
    #(12): "I1_Spitzer", "I2_Spitzer", "W2_WISE", "I3_Spitzer", "I4_Spitzer",
    #(17): "W3_WISE", "W4_WISE", "P70_Herschel", "P100_Herschel",
    #(21): "P160_Herschel", "S250_Herschel", "S350_Herschel", "S450_JCMT",
    #(25): "S500_Herschel", "S850_JCMT", "Band9_ALMA", "Band8_ALMA",
    #(29): "Band7_ALMA", "Band6_ALMA", "Band5_ALMA", "Band4_ALMA"

    writeon = True
    bands = range(0,11)
    #xyz distance
    for g in range(0,len(x_cen)):
        d_all = np.sqrt(pow(x - x_cen[g], 2.0) + pow(y - y_cen[g], 2.0) + pow(z - z_cen[g], 2.0))
        ind   = np.where((d_all < 2))
        massin = mstar_tot[ind]
        sfrin = sfr_tot[ind]
        xin = x[ind] - x_cen[g]
        yin = y[ind] - y_cen[g]
        zin = z[ind] - z_cen[g]
        vxin = vx[ind] - vx_cen[g]
        vyin = vy[ind] - vy_cen[g]
        vzin = vz[ind] - vz_cen[g]
        typegin = typeg[ind]
        id_halo_treein = id_halo_tree[ind]
        apssed = ap_mags[:,ind]
        abssed = ab_mags[:,ind]
        apssed = apssed[:,0,:]
        abssed = abssed[:,0,:]
       
        print(apssed.shape)
        if(writeon):
           f = open('cluster_z' + str(zlist[index]) + "_" + str(g) + '.txt', 'w')
           f.write("Mvir[Msun]: %s\n" % (str(mvir_cen[g])))
           f.write("#central galaxy has positions=0 and velocities=0\n")
           f.write("#all galaxies within a sphere of radius 2 comoving Mpc included\n")
           f.write("#mstar[Msun] sfr[Msun/yr] x[cMpc] y[cMpc] z[cMpc] vx[km/s] vy[km/s] vz[km/s] typeg(=0 centrals) id_halo_tree ap_mag[AB] (FUV_GALEX, NUV_GALEXi, u_SDSS, g_SDSS, r_SDSS, i_SDSS, z_SDSS, Y_VISTA, J_VISTA, H_VISTA, K_VISTA) ab_mag[AB] (FUV_GALEX, NUV_GALEXi, u_SDSS, g_SDSS, r_SDSS, i_SDSS, z_SDSS, Y_VISTA, J_VISTA, H_VISTA, K_VISTA)\n")
           for i in range(0,len(massin)):
               srt_to_write = str(massin[i]) + ' ' + str(sfrin[i]) + ' ' + str(xin[i]) + ' ' + str(yin[i]) +' ' + str(zin[i]) +' ' + str(vxin[i]) +' ' + str(vyin[i]) +' ' + str(vzin[i]) + ' ' +str(typegin[i]) + ' ' + str(id_halo_treein[i])
               for b in bands:
                   srt_to_write += ' ' + str(apssed[b,i])
               for b in bands:
                   srt_to_write += ' ' + str(abssed[b,i])
               srt_to_write += "\n"
               f.write(srt_to_write)
           f.close()
        


def main(model_dir, output_dir, redshift_table, subvols, obs_dir):


    plt = common.load_matplotlib()

    zlist = (5.0, 6.0)

    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge', 'mvir_hosthalo', 'sfr_disk',
                           'sfr_burst', 'position_x', 'position_y', 'position_z', 'vvir_hosthalo',
                           'velocity_x', 'velocity_y', 'velocity_z', 'id_halo_tree')}

    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5"
    fields_sed_ap = {'SED/ap_dust': ('total','disk'),}
    fields_sed_ab = {'SED/ab_dust': ('total','disk'),}

    for index, snapshot in enumerate(redshift_table[zlist]):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        seds_ap = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_ap, subvols, file_hdf5_sed)
        seds_ab = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_ab, subvols, file_hdf5_sed)

        prepare_data(hdf5_data, seds_ap, seds_ab, index, zlist)

if __name__ == '__main__':
    main(*common.parse_args())
