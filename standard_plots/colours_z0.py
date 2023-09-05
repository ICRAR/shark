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

def plot_colours(plt, outdir, obsdir, colours_dist):

    #plot colour distributions
    fig = plt.figure(figsize=(9.7,11.7))
    xtit = "$\\rm (u-r)$"
    ytit = "$\\rm dp/d(u-r)$"
    xmin, xmax, ymin, ymax = 0, 3.5, 0, 3
    xleg = xmax - 0.4 * (xmax - xmin)
    yleg = ymax - 0.15 * (ymax - ymin)
   
    subplots = (321, 322, 323, 324, 325, 326)
    indeces = (0, 1, 2, 3, 4, 5)
    labels  = ("(-17.13,-17.88)","(-17.88,-18.63)","(-18.63,-19.38)","(-19.38,-20.13)","(-20.13,-20.88)","(-20.88,-21.63)")
   
    file = obsdir+'/Colours/ur_colours_SDSS.data'
    cSDSS,dpSDSS = np.loadtxt(file,usecols=[0,1],unpack=True)
    columns = [0,38,77,112,151,190,229]
   
    for subplot, idx in zip(subplots, indeces):
   
        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, labels[idx])
        
        if(idx == 0):
           c1 = columns[idx]
        if(idx > 0):
           c1 = columns[idx]+1
        # Observed CDF
        yzeros = np.zeros(shape = len(cSDSS[c1:columns[idx+1]]))
        ax.plot(cSDSS[c1:columns[idx+1]], dpSDSS[c1:columns[idx+1]],'grey', label='SDSS')
        ax.fill_between(cSDSS[c1:columns[idx+1]], dpSDSS[c1:columns[idx+1]], yzeros, facecolor='grey', alpha=1,interpolate=True)
        
        # Predicted CDF
        y = colours_dist[idx,0,:]
        ind = np.where(y > 0.)
        ax.plot(xc[ind],y[ind],'r', label='Shark')
        
        common.prepare_legend(ax, ['grey','r'], loc = 2)

    common.savefig(outdir, fig, 'ur_colour_z0.pdf')

    fig = plt.figure(figsize=(9.7,11.7))
    xtit = "$\\rm (g-r)$"
    ytit = "$\\rm dp/d(g-r)$"
    xmin, xmax, ymin, ymax = 0, 1.1, 0, 6.5
    xleg = xmax - 0.4 * (xmax - xmin)
    yleg = ymax - 0.15 * (ymax - ymin)
    
    file = obsdir+'/Colours/gr_colours_SDSS.data'
    cSDSS,dpSDSS = np.loadtxt(file,usecols=[0,1],unpack=True)
    columns = [0,21,43,65,87,109,131]
    
    for subplot, idx in zip(subplots, indeces):
    
        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, labels[idx])
    
        if(idx == 0):
           c1 = columns[idx]
        if(idx > 0):
           c1 = columns[idx]+1
        # Observed CDF
        yzeros = np.zeros(shape = len(cSDSS[c1:columns[idx+1]]))
        ax.plot(cSDSS[c1:columns[idx+1]], dpSDSS[c1:columns[idx+1]],'grey', label='SDSS')
        ax.fill_between(cSDSS[c1:columns[idx+1]], dpSDSS[c1:columns[idx+1]], yzeros, facecolor='grey', alpha=1,interpolate=True)
    
        # Predicted CDF
        y = colours_dist[idx,1,:]
        ind = np.where(y > 0.)
        ax.plot(xc[ind],y[ind],'r', label='Shark')
    
        common.prepare_legend(ax, ['grey','r'], loc = 2)
    
    common.savefig(outdir, fig, 'gr_colour_z0.pdf')
  

def prepare_data(hdf5_data, phot_data, colours_dist, nbands):
   
    (h0, _, mdisk, mbulge, mhalo, mshalo, typeg, age, 
     sfr_disk, sfr_burst, id_gal) = hdf5_data
   
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
        for i in range(0,nbands):
            SEDs_dust[:,c,i] = indust[i,:]
        p = p + 1

    uband = 2
    gband = 3
    rband = 4
    ubandl = SEDs_dust[:,4,uband]
    gbandl = SEDs_dust[:,4,gband]
    rbandl = SEDs_dust[:,4,rband]

    for mag in range(0,len(magbins)-1):
        ind = np.where((ubandl < -1) & (ubandl > -30) & (gbandl < -1) & (gbandl > -30) & (rbandl < magbins[mag]) & (rbandl >= magbins[mag+1]))
        H, bins_edges  = np.histogram(ubandl[ind] - rbandl[ind],bins=np.append(cbins,cupp))
        colours_dist[mag,0,:] = colours_dist[mag,0,:] + H
        colours_dist[mag,0,:] = colours_dist[mag,0,:] / (len(ubandl[ind]) * dc)
        H, bins_edges  = np.histogram(gbandl[ind] - rbandl[ind],bins=np.append(cbins,cupp))
        colours_dist[mag,1,:] = colours_dist[mag,1,:] + H
        colours_dist[mag,1,:] = colours_dist[mag,1,:] / (len(gbandl[ind]) * dc)
 
def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()

    Variable_Ext = True
    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5"

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age',
                           'sfr_disk', 'sfr_burst', 'id_galaxy')}

    hdf5_data = common.read_data(model_dir, redshift_table[0], fields, subvols)

    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}

    if(Variable_Ext == False):
       seds = common.read_photometry_data(model_dir, redshift_table[0], fields_sed, subvols)
    else:
       seds = common.read_photometry_data_variable_tau_screen(model_dir, redshift_table[0], fields_sed, subvols, file_hdf5_sed)
 
    nbands = len(seds[0]) 

    colours_dist = np.zeros(shape = (len(magbins)-1, 2, len(cbins)))

    prepare_data(hdf5_data, seds, colours_dist, nbands)

    if(Variable_Ext):
       outdir = os.path.join(outdir, 'eagle-rr14')

    plot_colours(plt, outdir, obsdir, colours_dist)

if __name__ == '__main__':
    main(*common.parse_args())
