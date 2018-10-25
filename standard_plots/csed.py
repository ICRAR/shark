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

def plot_csed(plt, outdir, obsdir, h0, CSED, nbands):

    #wavelength in angstroms.
    file = obsdir+'/Models/Shark_SED_bands.dat'
    lambda_bands = np.loadtxt(file,usecols=[0],unpack=True)
    freq_bands   = c_light / (lambda_bands * 1e-10) #in Hz
    lambda_bands = np.log10(lambda_bands)

    print 'freq_bands', freq_bands
    xtit="$\\rm log_{10}(\lambda/Ang\, (rest-frame))$"
    ytit="$\\rm log_{10}(\\nu \\epsilon_{\\rm int}/ h\,W\, Mpc^{-3})$"

    fig = plt.figure(figsize=(7,12))

    subplots = (411, 412, 413, 414)
    idx = (0, 1, 2, 3)
    labels= ('z=0', 'z=0.25', 'z=0.5', 'z=1')
    obs = ('z0p05', 'z0p25', 'z0p5', 'z0p95')
    colors = ('Indigo','purple','Navy','MediumBlue','Green','MediumAquamarine','LightGreen','YellowGreen','Gold','Orange','Coral','OrangeRed','red','DarkRed','FireBrick','Crimson','IndianRed','LightCoral','Maroon','brown','Sienna','SaddleBrown','Chocolate','Peru','DarkGoldenrod','Goldenrod','SandyBrown')

    for subplot, idx in zip(subplots, idx):
        xmin, xmax, ymin, ymax = 3.0, 7.0, 30, 36
        xleg = xmin + 0.1 * (xmax-xmin)
        yleg = ymin + 0.1 * (ymax-ymin)

        ax = fig.add_subplot(subplot)
        if (idx == 3):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytit, locators=(1, 1, 1, 1))
        ax.text(xleg,yleg, labels[idx], fontsize=12)

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

        for xi,yi,c in zip(lambda_bands,np.log10(CSED[idx,4,:]*freq_bands),colors):
            ax.plot(xi,yi-np.log10(h0), 'x', markersize=6, color=c)
        ax.plot(lambda_bands,np.log10(CSED[idx,3,:]*freq_bands)-np.log10(h0), marker = 'o', mec = 'b', markersize=3, linewidth=1)
        ax.plot(lambda_bands,np.log10(CSED[idx,1,:]*freq_bands)-np.log10(h0), marker = 'd', mec = 'r', markersize=3, linewidth=1)
        ax.plot(lambda_bands,np.log10(CSED[idx,0,:]*freq_bands)-np.log10(h0), marker = 'p', mec = 'LightSalmon', markersize=3, linewidth=1)
        ax.plot(lambda_bands,np.log10(CSED[idx,4,:]*freq_bands)-np.log10(h0), 'k', linewidth=1)


        print lambda_bands,np.log10(CSED[idx,4,:]*freq_bands)-np.log10(h0)

    common.savefig(outdir, fig, "CSED_Shark.pdf")

def prepare_data(hdf5_data, phot_data, ids_sed, CSED, nbands, index):
   
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
    SEDs_dust   = phot_data[:,1,0,:,:]
    SEDs_nodust = phot_data[:,0,0,:,:]

    for i in range(0,nbands):
        for c in range(0,5):
            #calculate LF with bands with dust
            ind = np.where(SEDs_dust[:,c,i] < -1)
            #W Hz-1
            total = np.sum(pow(10.0,(SEDs_dust[ind,c,i]+48.6)/(-2.5))) * (4.0 * PI * pow(10.0*3.086e18,2.0)) * 1e-7 
            CSED[index,c,i] = total 

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age', 
                           'sfr_disk', 'sfr_burst', 'id_galaxy')}

    z = (0, 0.25, 0.5, 1) #, 1.0, 1.5, 2.0)
    snapshots = redshift_table[z]

    # Create histogram
    for index, snapshot in enumerate(snapshots):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        #sfh, delta_t, LBT = common.read_sfh(model_dir, snapshot, sfh_fields, subvols)
        seds, ids, nbands = common.read_photometry_data(model_dir, snapshot, subvols)
        
        if(index == 0):
            CSED = np.zeros(shape = (len(z), 5, nbands))

        prepare_data(hdf5_data, seds, ids, CSED, nbands, index)

        h0, volh = hdf5_data[0], hdf5_data[1]
        if(volh > 0.):
            CSED[index,:]   = CSED[index,:] / volh * pow(h0,3.0)

    # Take logs
    plot_csed(plt, outdir, obsdir, h0, CSED, nbands)
 
if __name__ == '__main__':
    main(*common.parse_args())
