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
"""SMF plots"""

import collections
import functools
import logging
import math

import numpy as np

import common
import utilities_statistics as us


logger = logging.getLogger(__name__)

##################################
# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
MpcToKpc = 1e3

##################################
# Mass function initialization
mlow = 4
mupp = 12
dm = 0.25
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

#velocity initialization
vlow = 1
vupp = 3.5
dv = 0.25
vbins = np.arange(vlow,vupp,dv)
xvf = vbins + dv/2.0

def plot_lf_z(plt, outdir, obsdir, hist_lf, co_vel_scaling):

    fig = plt.figure(figsize=(9.7,11.7))
    xtit = "$\\rm log_{10} (\\rm L_{\\rm CO}/\\rm Jy\\, km/s\\, Mpc^{-2})$"
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm L_{\\rm CO}}/{\\rm Mpc}^{-3})$"
    xmin, xmax, ymin, ymax = 4, 10, -6, -1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (321, 322, 323, 324, 325, 326)
    indeces = (0, 1, 2, 3, 4, 5)
    zs = (0, 0.5, 1, 1.5, 2, 3)

    colors   = ('purple','Navy','DarkTurquoise', 'Aquamarine', 'Green','Gold','Yellow','Orange','red','Chocolate')
    labels = ('1-0', '2-1','3-2','4-3','5-4','6-5','7-6','8-7','9-8','10-9') 

    for subplot, idx, z in zip(subplots, indeces, zs):

        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(z))

        # Predicted SMF
        for i in range(10):
            y = hist_lf[idx,i,:]
            ind = np.where(hist_lf[idx,i,:] != 0.)
            y = hist_lf[idx,i,ind]
            ax.plot(xmf[ind],y[0],color=colors[i], label=labels[i] if idx == 0 else None)
        if idx == 0:
           common.prepare_legend(ax, colors)

    common.savefig(outdir, fig, 'colf_z.pdf')


def prepare_data(hdf5_data, hdf5_data_gal, index, hist_lf):

    (h0, volh, mdisk, mbulge, rgas_disk, rgas_bulge, jgas_disk, jgas_bulge)  = hdf5_data_gal
    (co_disk, co_bulge) = hdf5_data


    vdisk = jgas_disk / rgas_disk / 2.0  #in km/s
    vbulge = jgas_bulge / rgas_bulge / 2.0 #in km/s

    co_total = co_disk + co_bulge
    print co_total.shape, vdisk.shape
     
    ind = np.where( (mdisk +  mbulge) > 0) 
    v_equiv  = (co_disk[:,0] * vdisk[ind] + co_bulge[:,0] * vbulge[ind]) / (co_disk[:,0] + co_bulge[:,0])
    print v_equiv, co_disk[:,0] + co_bulge[:,0]
    print 'list of gals at index', index
    if(index  == 4):
       for a,b in zip(co_total[:,0], v_equiv[:]):
           if(a > 0 and a < 1e15 and b > 0 and b < 1e5):
              print a,b

    for i in range(10):
        ind = np.where(co_total[:,i] > 0.0)
        H, _ = np.histogram(np.log10(co_total[ind,i]),bins=np.append(mbins,mupp))
        hist_lf[index,i,:] = hist_lf[index,i,:] + H
        hist_lf[index,i,:] = hist_lf[index,i,:]/volh/dm

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist = (0, 0.5, 1, 1.5, 2, 3)

    plt = common.load_matplotlib()

    # Histograms
    hist_lf        = np.zeros(shape = (len(zlist), 10, len(mbins)))
    co_vel_scaling = np.zeros(shape = (len(zlist), 2, len(vbins)))

    fields = {'galaxies': ('LCO_disk', 'LCO_bulge')}
    fields_gal = {'galaxies': ('mstars_disk', 'mstars_bulge','rgas_disk', 'rgas_bulge', 'specific_angular_momentum_disk_gas', 'specific_angular_momentum_bulge_gas')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_co_data(modeldir, snapshot, fields, subvols)
        hdf5_data_gal = common.read_data(modeldir, snapshot, fields_gal, subvols)
        prepare_data(hdf5_data, hdf5_data_gal, index, hist_lf, co_vel_scaling)


    # Take logs
    ind = np.where(hist_lf > 0.)
    hist_lf[ind] = np.log10(hist_lf[ind]) 

    write = True
    if write:
       for index, z in enumerate(zlist):
           with open('Shark-Lagos19-COLFs_%s.txt' % str(z), 'wb') as fil:
                fil.write("#Galaxies from Shark (Lagos et al. 2018, 2019) using CO modelling of Lagos et al. (2012)\n")
                fil.write("#CO LFs for lines (1-0) to (10-9)\n")
                fil.write("#Units of CO luminosity in [Jy km/s Mpc^2] and presented in log10\n")
                fil.write("#Units number density [Mpc^-3 dex^-1] and presented in log10\n")
                fil.write("#log10(LCO) logn1-0 logn2-1 logn3-2 logn4-3 logn5-4 logn6-5 logn7-6 logn8-7 logn9-8 logn10-9\n") 
                for a,b,c,d,e,f,g,h,i,j,k in zip(xmf[:],hist_lf[index,0,:], hist_lf[index,1,:],hist_lf[index,2,:],hist_lf[index,3,:],hist_lf[index,4,:],hist_lf[index,5,:],hist_lf[index,6,:],hist_lf[index,7,:],hist_lf[index,8,:],hist_lf[index,9,:]):
                    fil.write("%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n" % (a,b,c,d,e,f,g,h,i,j,k))


    plot_lf_z(plt, outdir, obsdir, hist_lf)

if __name__ == '__main__':
    main(*common.parse_args())
