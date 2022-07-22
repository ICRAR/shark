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


observation = collections.namedtuple('observation', 'label x y yerrup yerrdn err_absolute')

logger = logging.getLogger(__name__)

##################################
# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
MpcToKpc = 1e3

mlow = 7.0
mupp = 12.0
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

def prepare_data(hdf5_data, index, zlist, mzr, mzr_sf, mszr):

    (h0, volh, id_galaxy, sfr_disk, sfr_burst, mdisk, mbulge, rstar_disk, mBH, mHI, mH2, 
     mgas_disk, mHI_bulge, mH2_bulge, mgas_bulge, mgas_metals_disk, mgas_metals_bulge, 
     mstars_metals_disk, mstars_metals_bulge, typeg, mvir_hosthalo, rstar_bulge, 
     mbulge_mergers, mbulge_diskins, mbulge_mergers_assembly, mbulge_diskins_assembly) = hdf5_data


    bin_it = functools.partial(us.wmedians, xbins=xmf)

    zstar = (mstars_metals_disk + mstars_metals_bulge) / (mdisk + mbulge)
    rcomb = (rstar_disk * mdisk + rstar_bulge * mbulge) / (mdisk + mbulge) / h0 * 1e3
    mgas = mgas_disk+mgas_bulge
    mgas_metals = mgas_metals_disk + mgas_metals_bulge
    mass = np.zeros(shape = len(sfr_disk))
    ssfr = np.zeros(shape = len(sfr_disk)) 
    sfr = np.zeros(shape = len(sfr_disk))

    ind = np.where(mdisk + mbulge > 0)
    mass[ind] = np.log10((mdisk[ind] + mbulge[ind])/h0)
    ssfr[ind] = (sfr_disk[ind] + sfr_burst[ind]) / (mdisk[ind] + mbulge[ind]) #in Gyr^-1
    sfr[ind] = np.log10((sfr_disk[ind] + sfr_burst[ind])/h0/1e9)
   
    # calculate main sequence:
    ms = np.zeros(shape = len(xmf))
    for j in range(0,len(xmf)):
        ind = np.where((mass > xmf[j]-dm/2.0) & (mass <= xmf[j]+dm/2.0) & (sfr - mass + 9 > -5 + 0.5*zlist[index]) & (typeg == 0))
        if(len(mass[ind] > 0)):
               ms[j] = np.median(sfr[ind] - mass[ind] + 9.0)
    ind = np.where((ms != 0) & (xmf > 8.5) & (xmf < 9.5))
    (ms_fit_slope, ms_fit_offs) = np.polyfit(xmf[ind],ms[ind],1)

    dist_ms = (sfr - mass + 9.0) - (ms_fit_slope * mass + ms_fit_offs)
  
    ind = np.where((mass > 8.5) & (mass < 9.5) & (typeg == 0))
    print(np.median(dist_ms[ind]))

    ind = np.where((mgas_metals > 0.0) & (mgas > 0))
    mzr[index,:] = bin_it(x=mass[ind], y=np.log10((mgas_metals[ind]/mgas[ind]/Zsun)))

    ind = np.where((mgas_metals > 0.0) & (mgas > 0) & (dist_ms > -1) & (dist_ms < 1)) 
    mzr_sf[index,:] = bin_it(x=mass[ind], y=np.log10((mgas_metals[ind]/mgas[ind]/Zsun)))

    ind = np.where(mstars_metals_disk+mstars_metals_bulge > 0.0)
    mszr[index,:] = bin_it(x=mass[ind], y=np.log10(((mstars_metals_disk[ind]+mstars_metals_bulge[ind])/(mdisk[ind]+mbulge[ind])/Zsun)))

    #ind = np.where((mdisk+mbulge)/h0 > 1e8)
    #print('#id_galaxy sfr mstellar zstar r50 type_galaxy at z=%s' % str(zlist[index]))
    #for idg,a,b,c,d,e in zip(id_galaxy[ind], (sfr_disk[ind]+sfr_burst[ind])/h0/1e9, (mdisk[ind]+mbulge[ind])/h0, zstar[ind], rcomb[ind], typeg[ind]):
    #    print (idg,a,b,c,d,e,zlist[index])

    print('#Ms Zg delta_Zg_16 delta_Zg_84 Zg_sf delta_Zg_sf_16 delta_Zg_sf_84 Zs delta_Zs_16 delta_Zs_84 at z=%s' % str(zlist[index]))
    for a,b,c,d,e,f,g,h,i,j in zip(xmf[:], mzr[index,0,:], mzr[index,1,:], mzr[index,2,:], mzr_sf[index,0,:], mzr_sf[index,1,:], mzr_sf[index,2,:], mszr[index,0,:], mszr[index,1,:], mszr[index,2,:]):
        print(a,b,c,d,e,f,g,h,i,j)

    return (h0, volh)
def plot_mzr_z0(plt, outdir, obsdir, h0, mzr, mzr_sf, mszr, zlist):

    fig = plt.figure(figsize=(4.5,10))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\\rm Z_{\\rm gas}/Z_{\odot})$"
    xmin, xmax, ymin, ymax = 8, 12, -2, 1

    ax = fig.add_subplot(311)
    plt.subplots_adjust(bottom=0.15, left=0.18)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, ' ', ytit, locators=(0.1, 1, 0.1))

    #MZR z=0
    corrzsun = 8.69 #solar oxygen abundance in units of 12 + log(O/H)
    hobs = 0.72
    #add cosmology correction plus IMF correction that goes into the stellar mass.
    corr_cos = np.log10(pow(hobs,2)/pow(h0,2)) - 0.09

    def load_obs_as(ax):
         lm, mz, mzdn, mzup = common.load_observation(obsdir, 'MZR/MMAdrews13.dat', [0,1,2,3])
         hobs = 0.7
         #add cosmology correction plus IMF correction that goes into the stellar mass.
         corr_cos = np.log10(pow(hobs,2)/pow(h0,2)) - 0.09
         common.errorbars(ax, lm+ corr_cos, mz - corrzsun, mzdn - corrzsun, mzup - corrzsun, 'grey', 's', label='Andrews+13')
         #correction for Tremonti is the same.
         lm, mz, mzdn, mzup = common.load_observation(obsdir, 'MZR/Tremonti04.dat', [0,1,2,3])
         common.errorbars(ax, lm+ corr_cos, mz - corrzsun, mzdn - corrzsun, mzup - corrzsun, 'grey', 'o', label="Tremonti+04")

    load_obs_as(ax)
    cols = ['DarkRed','Red','Crimson','OrangeRed','DarkOrange','Gold','Yellow','Chocolate','GoldenRod','SandyBrown','Green','YellowGreen','LawnGreen','Aqua','DarkTurquoise','Teal','DeepSkyBlue','Blue','Navy']

    def load_predictions(ax, zlist, mzr, cols):
         for i, z in enumerate(zlist):
             ind = np.where(mzr[i,0,:] != 0)
             yplot = (mzr[i,0,ind])
             errdn = (mzr[i,1,ind])
             errup = (mzr[i,2,ind])
             xplot = xmf[ind]
             if(z == 0 or z==1 or z == 3):
                ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=cols[i], alpha=0.4,interpolate=True)
                ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=cols[i], alpha=0.4,interpolate=True)
             ax.plot(xplot,yplot[0],color=cols[i], linestyle='solid') #, label = 'z=%s' % str(z))
   
    load_predictions(ax,zlist,mzr, cols) 
    common.prepare_legend(ax, ['grey','grey','grey'], loc=4)

    ax = fig.add_subplot(312)
    plt.subplots_adjust(bottom=0.15, left=0.18)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, ' ', ytit, locators=(0.1, 1, 0.1))

    load_obs_as(ax)
    load_predictions(ax,zlist,mzr_sf,cols)

    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\\rm Z_{\\star}/Z_{\odot})$"
    xmin, xmax, ymin, ymax = 8, 12, -2, 1

    ax = fig.add_subplot(313)
    plt.subplots_adjust(bottom=0.15, left=0.18)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    #MZR z=0
    lm, mz, mzdn, mzup = common.load_observation(obsdir, 'MZR/MSZR-Gallazzi05.dat', [0,1,2,3])
    common.errorbars(ax, lm[0:7], mz[0:7], mzdn[0:7], mzup[0:7], 'grey', 'D', label='Kirby+13')
    common.errorbars(ax, lm[7:22], mz[7:22], mzdn[7:22], mzup[7:22], 'grey', 'o', label='Gallazzi+05')

    load_predictions(ax,zlist,mszr,cols)

    common.prepare_legend(ax, ['grey','grey'], loc=4)

    common.savefig(outdir, fig, 'mzr_evolution.pdf')


def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist = (0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)

    plt = common.load_matplotlib()

    mzr         = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mzr_sf     = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mszr        = np.zeros(shape = (len(zlist), 3, len(xmf)))

    fields = {'galaxies': ('id_galaxy', 'sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge',
                           'mgas_metals_disk', 'mgas_metals_bulge',
                           'mstars_metals_disk', 'mstars_metals_bulge', 'type', 
			   'mvir_hosthalo', 'rstar_bulge', 'mstars_burst_mergers', 
                           'mstars_burst_diskinstabilities', 'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        (h0, volh) = prepare_data(hdf5_data, index, zlist, mzr, mzr_sf, mszr)

    plot_mzr_z0(plt, outdir, obsdir, h0, mzr, mzr_sf, mszr, zlist)

if __name__ == '__main__':
    main(*common.parse_args())
