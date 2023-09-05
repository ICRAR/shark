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

import common


##################################
mlow = 10
mupp = 15
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

# Constants
GyrToYr = 1e9
zsun = 0.0189
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3
c_light = 299792458.0 #m/s

def plot_fraction_hydro(plt, outdir, obsdir, f_q, z):

    print("Will make plot now")
    xtit="$\\rm log_{10}(M_{\\rm halo}/M_{\odot})$"
    ytit="$\\rm f_{\\rm hot-halo}$"

    xmin, xmax, ymin, ymax = 10, 15.0, -0.05, 1.05
    xleg = xmin + 0.1 * (xmax-xmin)
    yleg = ymax - 0.07 * (ymax-ymin)

    fig = plt.figure(figsize=(5,5))
    colors = ('Navy','Blue','RoyalBlue','SkyBlue','LightSalmon','IndianRed','Crimson','Red','DarkRed')

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(2, 2, 0.15, 0.15))
    #ax.text(xleg,yleg, labelages[idx])

    for i,j in enumerate(z):
        ind = np.where(f_q[i,:] !=0)
        xin = xmf[ind]
        yin = f_q[i,ind]
        ax.plot(xin, yin[0], color=colors[i], linewidth=3, label = 'z=%s' % str(j))

    common.prepare_legend(ax, ['k', 'k'], loc=2)
    plt.tight_layout()
    common.savefig(outdir, fig, "hothalo_fraction_z.pdf")
    

def prepare_data(hdf5_data, hdf5_data_halo, index, f_q, read_hydroeq):
   
    if(read_hydroeq):
       (h0, _, mdisk, mbulge, mhalo, mshalo, typeg, age, 
        sfr_disk, sfr_burst, id_gal, mbh, mhot, mreheated,
        on_hydrostatic_eq, id_halo) = hdf5_data
       (h0, _, halo_id, age_80, age_50) = hdf5_data_halo
    else:
        (h0, _, mdisk, mbulge, mhalo, mshalo, typeg, age,
        sfr_disk, sfr_burst, id_gal, mbh, mhot, mreheated) = hdf5_data

    ind = np.where((typeg ==0) & (mhalo/h0 > 10**13.585662))
    id_halosin = id_halo[ind]
    ages_h = np.zeros ( shape = (2, len(typeg[ind])))
    j = 0
    i = 0

    while (j < len(typeg[ind])):
        if(halo_id[i] == id_halosin[j]):
           ages_h[0,j] = age_80[i]
           ages_h[1,j] = age_50[i]
           j = j + 1
        i = i + 1

    print(np.median(ages_h[0,:]),  np.median(ages_h[1,:]))

    ind = np.where(((mdisk + mbulge)/h0 > 1e12) & (typeg ==0))
    sfrin = (sfr_disk[ind] + sfr_burst[ind])/h0/1e9
    sfr_fracin = sfr_burst[ind] / (sfr_disk[ind] + sfr_burst[ind])
    mhaloin = np.log10(mhalo[ind]/h0)
    typein = typeg[ind]
    id_halosin = id_halo[ind]
    mbhin = np.log10(mbh[ind]/h0)
    ages = np.zeros ( shape = (2, len(typein)))
    j = 0
    i = 0

    while (j < len(typein)): 
        if(halo_id[i] == id_halosin[j]):
           ages[0,j] = age_80[i]
           ages[1,j] = age_50[i]
           j = j + 1
        i = i + 1 

    ind = np.where(sfrin > 20)
    print(np.median(mhaloin[ind]), np.median(typein[ind]), np.median(ages[0,ind]),  np.median(ages[1,ind]), np.median(mbhin[ind]), np.median(sfr_fracin[ind]))
    ind = np.where(sfrin < 20)
    print(np.median(mhaloin[ind]), np.median(typein[ind]), np.median(ages[0,ind]),  np.median(ages[1,ind]), np.median(mbhin[ind]), np.median(sfr_fracin[ind]))


    if(read_hydroeq):
       #ind = np.where((typeg == 0) & (mhalo/h0 > 3e12))
       #on_hydrostatic_eq[ind] = 1
       for i,m in enumerate(xmf):
           ind = np.where((typeg == 0) & (on_hydrostatic_eq >= 0) & (mhalo/h0 > 10**(m-dm/2.0)) & (mhalo/h0 < 10**(m+dm/2.0)))
           n_all = len(typeg[ind])
           if(n_all > 9):
              ind = np.where((typeg == 0) & (on_hydrostatic_eq == 1) & (mhalo/h0 > 10**(m-dm/2.0)) & (mhalo/h0 < 10**(m+dm/2.0)))
              n_hydro = len(typeg[ind])
              f_q[index,i] = (n_hydro + 0.0) / (n_all + 0.0)
           #print(m, n_hydro, n_all)
   
   
       ind =np.where((typeg == 0) & (mhalo > 1e14))
       #print("massive halos", (on_hydrostatic_eq[ind]))

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    read_hydroeq = True

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()
    if(read_hydroeq):
       fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                              'mvir_subhalo', 'type', 'mean_stellar_age', 
                              'sfr_disk', 'sfr_burst', 'id_galaxy', 'm_bh',
                              'mhot','mreheated','on_hydrostatic_eq', 'id_halo_tree')}
       fields_halo = {'halo': ('halo_id', 'age_80', 'age_50')}
    else:
        fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                              'mvir_subhalo', 'type', 'mean_stellar_age',
                              'sfr_disk', 'sfr_burst', 'id_galaxy', 'm_bh',
                              'mhot','mreheated')}


    z = [0, 0.5, 1, 2, 3, 4, 6]
    snapshots = redshift_table[z]

    f_q = np.zeros(shape = (len(z), len(xmf)))

    # Create histogram
    for index, snapshot in enumerate(snapshots):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        hdf5_data_halo = common.read_data(model_dir, snapshot, fields_halo, subvols)

        prepare_data(hdf5_data, hdf5_data_halo, index, f_q, read_hydroeq)

    if(read_hydroeq):
       plot_fraction_hydro(plt, outdir, obsdir, f_q, z)


  
if __name__ == '__main__':
    main(*common.parse_args())
