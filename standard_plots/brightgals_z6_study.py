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

def plot_distribution(plt, outdir, x, y, mag):

    fig = plt.figure(figsize=(5,4.5))

    # Total ##################################
    xtit="x/Mpc"
    ytit="y/Mpc"
    xmin, xmax, ymin, ymax = -0.7, 0.7, -0.7, 0.7

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.2)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 0.2, 0.1))

    inp = np.where((mag < -10) & (mag > -30))
    print(min(abs(mag[inp]/15.0)),max(abs(mag[inp]/15.0)))
    xin = x[inp]
    yin = y[inp]
    sizes = abs(mag[inp]/10.0)
    mags = mag[inp]

    sortedm = np.argsort(abs(mags))
    mags = mags[sortedm]
    sizes = sizes[sortedm]
    xin = xin[sortedm]
    yin = yin[sortedm]

    colors = ['PowderBlue', 'LightSkyBlue','DeepSkyBlue', 'Blue', 'DarkBlue']
    print(max(mag[inp]))
    magscuts = [-15,-17.5,-18.5,-20,-30]
    for i in range(0,len(xin)):
        marker = 'o'
        markersize = sizes[i]
        if(mags[i] > -15):
           col = colors[0]
        if((mags[i] <= -15) & (mags[i] >= -17.5)):
           col = colors[1]
        if((mags[i] <= -17.5) & (mags[i] >= -18.5)):
           col = colors[2]
        if((mags[i] <= -18.5) & (mags[i] >= -19.8)):
           col = colors[3]
        if(mags[i] <= -19.8):
           print('bright gal')
           col = colors[4]
           marker='*'
           markersize = 8

        ax.plot(xin[i],yin[i], color=col, marker=marker, markersize = markersize, linewidth=0)

    xs =[-0.5,0.5,0.5,-0.5,-0.5]
    ys = [-0.5,-0.5,0.5,0.5,-0.5]
    ax.plot(xs,ys,linestyle='solid',color='r')
    xs =[-0.1,0.1,0.1,-0.1,-0.1]
    ys = [-0.1,-0.1,0.1,0.1,-0.1]
    ax.plot(xs,ys,linestyle='dotted',color='r')

    common.savefig(outdir, fig, 'clustering_around_protoclusterz6.pdf')

def prepare_data(hdf5_data, seds_ab, hdf5_data_halo, index, zlist):


    # Unpack data
    (h0, _, typeg, mdisk, mbulge, _, _, mHI, mH2, mgas,
     mHI_bulge, mH2_bulge, mgas_bulge, mvir, sfrd, sfrb, 
     x, y, z, vvir, vx, vy, vz, id_halo_tree) = hdf5_data
    (h0, _, mvirhaloz0, mvirhalo) = hdf5_data_halo
    ids_halocat =range(0,len(mvirhalo))

    mvirhaloz0_galcat = np.zeros(shape = (len(typeg)))
    ind = np.where(typeg == 0)
    g = 0
    for i in range(0,len(x)):
        if(typeg[i] == 0):
           mvirhaloz0_galcat[i] = mvirhaloz0[g]
           g = g + 1

    ab_mags = seds_ab[0]

    uvmag = ab_mags[0]
    ind = np.where( (mdisk +  mbulge) > 0)
    mvir = mvir[ind]
    vvir = vvir[ind]
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
    mvirhaloz0_galcat = mvirhaloz0_galcat[ind]
    
    #mm = np.where(mvir == max(mvir))
    #print("maximum mass", mvir[mm], mvirhaloz0_galcat[mm], max(mvirhaloz0_galcat))

    d = 0.2 * (1 + 6.0) #comoving Mpc
    rvir  = G * mvir / pow(vvir,2.0) / h0

    #select L* galaxies
    lstar = np.where((uvmag < -19.8) & (uvmag > -30))
    xlstar = x[lstar]
    ylstar = y[lstar]
    zlstar = z[lstar]

    mvirlstar = mvir[lstar]
    rvirlstar = rvir[lstar]
    id_halo_tree_lstar = id_halo_tree[lstar]
    typeg_lstar = typeg[lstar]
    mvirz0_lstar = mvirhaloz0_galcat[lstar]

    #xselec, yselect, zselect for the galaxy we'll select
    distance_muse = 0.5 * (1 + 6.0) * h0

    mvir_highcluster = []
    rvir_highcluster = []
    n_neigh_highcluster = []
    ids_highcluster = np.array([])
    mvirz0_highcluster = []
    for i in range(0,len(xlstar)):
        d_bright = np.sqrt(pow(xlstar - xlstar[i], 2.0) + pow(ylstar - ylstar[i], 2.0))
        n_neigh = np.count_nonzero((d_bright <= d) & (abs(zlstar-zlstar[i]) < 2))
        neigh = np.where(d_bright <= d)
        if(n_neigh >= 3):
           uniqueids = True
           if(len(ids_highcluster) > 0):
              allids = np.append(ids_highcluster, id_halo_tree_lstar[neigh])
           else:
              allids = id_halo_tree_lstar[neigh]
           if(len(ids_highcluster) > 0):
              nonrep = np.unique(allids)
              #if all IDs are repeated, then don't count this as anew overdensity of interest
              if(len(nonrep) <= len(allids)-n_neigh):
                 uniqueids = False
           if(uniqueids): 
              mvir_highcluster.append(mvirlstar[i])
              rvir_highcluster.append(rvirlstar[i])
              mvirz0_highcluster.append(mvirz0_lstar[i])
              n_neigh_highcluster.append(n_neigh)
              ids_highcluster = allids

              #select bright galaxies within muse pointing
              box = np.where((abs(xlstar-xlstar[i]) < distance_muse) & (abs(ylstar-ylstar[i]) < distance_muse) & ( abs(zlstar-zlstar[i]) < distance_muse*2.0))
              neigh_largebox = len(xlstar[box])
              print("Mvir:", mvirlstar[i]/h0, n_neigh, id_halo_tree_lstar[i],mvirz0_lstar[i]/h0, neigh_largebox)
              if(mvirz0_lstar[i]/h0 > 9.24e14):
                xselec = xlstar[i]
                yselec = ylstar[i]
                zselec = zlstar[i]

    cen = np.where(typeg == 0)
    mvir_mostmass = 1.0/(np.sort(1.0/mvir[cen]))
    print(mvir_mostmass[9],mvir_mostmass[49],mvir_mostmass[99],mvir_mostmass[999],len(mvir_mostmass))
    print('median halo mass, highly clustered L*', np.mean(mvir_highcluster/h0), np.mean(rvir_highcluster),len(mvir_highcluster)) 

    distance_muse = 0.5 * (1 + 6.0) * h0
    box = np.where((abs(x-xselec) < distance_muse + 1) & (abs(y-yselec) < distance_muse + 1) & ( abs(z-zselec) < distance_muse))
    xneigh = (x[box]-xselec) / h0 / (1+6.0)
    yneigh = (y[box]-yselec) / h0 / (1+6.0)
    zneigh = (z[box]-zselec) / h0 / (1+6.0)
    uvmag = uvmag[box]
   
    return(xneigh, yneigh, uvmag)

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):


    plt = common.load_matplotlib()

    zlist = [6.0]

    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge', 'mvir_hosthalo', 'sfr_disk',
                           'sfr_burst', 'position_x', 'position_y', 'position_z', 'vvir_hosthalo',
                           'velocity_x', 'velocity_y', 'velocity_z', 'id_halo_tree')}

    fields_halo = {'halo': ('final_z0_mvir','mvir')}
    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5"
    fields_sed_ab = {'SED/ab_dust': ('total','disk'),}

    for index, snapshot in enumerate(redshift_table[zlist]):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        hdf5_data_halo = common.read_data(model_dir, snapshot, fields_halo, subvols)

        seds_ab = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_ab, subvols, file_hdf5_sed)

        (x, y, mag) = prepare_data(hdf5_data, seds_ab, hdf5_data_halo, index, zlist)

    plot_distribution(plt, output_dir, x, y, mag)

if __name__ == '__main__':
    main(*common.parse_args())
