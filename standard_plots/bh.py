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
"""Size plots"""

import functools

import numpy as np

import os
import common
import utilities_statistics as us

# Initialize arguments
zlist = (0,  0.5,  1, 1.5,  2, 3, 4, 6)

##################################
#Constants
RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun

mlow = 4.5
mupp = 12.5
dm = 0.25
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

mdotlow = -5.0
mdotupp = 1.0
dmdot = 0.2
mdotbins = np.arange(mdotlow,mdotupp,dmdot)
xmdotf = mdotbins + dmdot/2.0


def prepare_data(hdf5_data, index, spinbh, spinmdot, mdotmbh, mdotmbh_hh, BHSFR, read_spin):

    if(read_spin == True):
       (h0, _, spin, mbh, mdot_hh, mdot_sb, sfr_disk, sfr_burst, mdisk, mbulge, typeg) = hdf5_data
       ind = np.where((np.isnan(spin) == True) | (np.isinf(spin) == True))
       if (len(spin[ind]) > 0):
         print("Number of galaxies with a spin of NaN:", len(spin[ind]))
       else:
         print("All galaxies have well defined BH spins")

    else:
       (h0, _, mbh, mdot_hh, mdot_sb, sfr_disk, sfr_burst, mdisk, mbulge, typeg) = hdf5_data

    print(np.log10(max(mbh)/h0))
    ind = np.where(mbh/h0 > 1e8)
    print("Number of galaxies with BH mass > 1e8Msun:", len(mbh[ind]), " at redshift", zlist[index])

    bin_it   = functools.partial(us.wmedians, xbins=xmf, nmin=5)
    bin_it_mdot = functools.partial(us.wmedians, xbins=xmdotf)

    if(read_spin == True):
       ind = np.where((np.isnan(spin) == True) | (np.isinf(spin) == True))
       if (len(spin[ind]) > 0):
            print("Number of galaxies with a BH spin of NaN:", len(spin[ind]))
       else:
            print("All galaxies have well defined BH spin")
   
       print(max(spin), min(spin))
       spin = np.abs(spin)
       ind = np.where(mbh > 0)
       spinbh[index,:] = bin_it(x=np.log10(mbh[ind]) - np.log10(float(h0)),
                               y=spin[ind])
       print(spinbh[index,:])
       ind = np.where((mbh > 0) & ((mdot_hh + mdot_sb) > 0 ))
       spinmdot[index,:] = bin_it_mdot(x=np.log10(mdot_hh[ind] + mdot_sb[ind]) - np.log10(float(h0) * 1e9),
                               y=spin[ind])
    ind = np.where((mbh > 0) & ((mdot_hh + mdot_sb) > 0 ))
    mdotmbh[index,:] = bin_it(x=np.log10(mbh[ind]) - np.log10(float(h0)),
                         y=np.log10(mdot_hh[ind] + mdot_sb[ind]) - np.log10(float(h0) * 1e9))

    mdotmbh_hh[index,:] = bin_it(x=np.log10(mbh[ind]) - np.log10(float(h0)),
                         y=np.log10(mdot_hh[ind]) - np.log10(float(h0) * 1e9))


    ssfr = (sfr_disk + sfr_burst) / 1e9 / (mbulge + mdisk)
    #apply lower limit to ssfr
    ind = np.where(ssfr < 1e-14)
    ssfr[ind] = 2e-14

    ind = np.where((mbh > 0) & (ssfr > 1e-14) & ((mbulge + mdisk)/h0 > 1e9) & (typeg <= 0))
    BHSFR[index,:] = bin_it(x=np.log10(mbh[ind]) - np.log10(float(h0)),
                            y=np.log10(ssfr[ind]))


def plot_spine_BH(plt, outdir, obsdir, spinbh, spinmdot, mdotmbh):

    fig = plt.figure(figsize=(5,9))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm BH}/M_{\odot})$"
    ytit = "spin"

    xmin, xmax, ymin, ymax = 4, 11, -0.05, 1.05
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(311)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'z=0')

    #spin bh mass relation
    ind = np.where(spinbh[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = spinbh[0,0,ind]
        errdn = spinbh[0,1,ind]
        errup = spinbh[0,2,ind]
        ax.plot(xplot,yplot[0],color='k',label="Shark")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='grey', interpolate=True)


    #mdot bh mass relation
    ytit = "$\\rm log_{10} (\\rm \\dot{M}_{\\rm BH}/M_{\\odot}\\, yr^{-1})$"

    ax = fig.add_subplot(312)
    plt.subplots_adjust(bottom=0.15, left=0.15)
    ymin, ymax = -4, 1

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'z=0')

    ind = np.where(mdotmbh[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = mdotmbh[0,0,ind]
        errdn = mdotmbh[0,1,ind]
        errup = mdotmbh[0,2,ind]
        ax.plot(xplot,yplot[0],color='k',label="Shark")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='grey', interpolate=True)

    common.prepare_legend(ax, ['k'], loc=2)


    #spin mdot relation
    xtit = "$\\rm log_{10} (\\rm \\dot{M}_{\\rm BH}/M_{\\odot}\\, yr^{-1})$"
    ytit = "spin"

    xmin, xmax, ymin, ymax = -4, 1, -0.05, 1.05
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(313)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'z=0')

    ind = np.where(spinmdot[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmdotf[ind]
        yplot = spinmdot[0,0,ind]
        errdn = spinmdot[0,1,ind]
        errup = spinmdot[0,2,ind]
        ax.plot(xplot,yplot[0],color='k',label="Shark")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='grey', interpolate=True)

    common.prepare_legend(ax, ['k'], loc=2)

    plt.tight_layout()
    common.savefig(outdir, fig, 'spin-BH.pdf')

def plot_macc_BH(plt, outdir, obsdir, mdotmbh, mdotmbh_hh, BHSFR):

    #cols=('red','LightSalmon','LimeGreen','DarkGreen','DarkTurquoise','blue')
    cols= ('Indigo','purple','Navy','DarkTurquoise', 'Aquamarine', 'Green','PaleGreen','GreenYellow','Gold','Yellow','Orange','OrangeRed','red') 
#,'DarkRed','FireBrick','Crimson','IndianRed','LightCoral','Maroon','brown','Sienna','SaddleBrown','Chocolate','Peru','DarkGoldenrod','Goldenrod','SandyBrown')


    fig = plt.figure(figsize=(5,7))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm BH}/M_{\odot})$"
    #mdot bh mass relation
    ytit = "$\\rm log_{10} (\\rm \\dot{M}_{\\rm BH,all}/M_{\\odot}\\, yr^{-1})$"

    ax = fig.add_subplot(211)
    plt.subplots_adjust(bottom=0.15, left=0.15)
    xmin, xmax, ymin, ymax = 4, 11, -5, 1

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    for i,c in enumerate(zlist):
        ind = np.where(mdotmbh[i,0,:] != 0)
        if(len(xmf[ind]) > 0):
            xplot = xmf[ind]
            yplot = mdotmbh[i,0,ind]
            errdn = mdotmbh[i,1,ind]
            errup = mdotmbh[i,2,ind]
            ax.plot(xplot,yplot[0],color=cols[i],label="z=%s" % str(c))
            if(i == 0 or i == 5):
               ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=cols[i], alpha=0.5, interpolate=True)
               ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=cols[i], alpha=0.5, interpolate=True)
    
    common.prepare_legend(ax, cols, loc=2)


    #spin mdot relation
    ytit = "$\\rm log_{10} (\\rm \\dot{M}_{\\rm BH,hh}/M_{\\odot}\\, yr^{-1})$"

    ax = fig.add_subplot(212)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    for i,c in enumerate(zlist):
        ind = np.where(mdotmbh_hh[i,0,:] != 0)
        if(len(xmf[ind]) > 0):
            xplot = xmf[ind]
            yplot = mdotmbh_hh[i,0,ind]
            errdn = mdotmbh_hh[i,1,ind]
            errup = mdotmbh_hh[i,2,ind]
            ax.plot(xplot,yplot[0],color=cols[i],label="z=%s" % str(c))
            if(i == 0 or i == 5 or i == 7):
               ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=cols[i], alpha=0.5, interpolate=True)
               ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=cols[i], alpha=0.5, interpolate=True)
    
    common.prepare_legend(ax, cols, loc=2)

    plt.tight_layout()
    common.savefig(outdir, fig, 'macc-BH.pdf')


    cols= ('Indigo','Navy','DarkTurquoise', 'Aquamarine', 'PaleGreen', 'Gold','Orange','red') 

    #SSFR vs BH mass
    fig = plt.figure(figsize=(3.5,4.5))
    ytit = "$\\rm log_{10} (\\rm sSFR/M_{\odot} yr^{-1})$"
    xtit = "$\\rm log_{10} (\\rm M_{\\rm BH}/M_{\odot})$"

    xmin, xmax, ymin, ymax = 5, 11, -14, -7
    xleg = xmax - 0.4 * (xmax - xmin)
    yleg = ymax + 0.05 * (ymax - ymin)

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(7, -6.9, 'Shark v2.0', fontsize=12)

    #Predicted BH-bulge mass relation
    for i,c in enumerate(zlist):
       ind = np.where(BHSFR[i,0,:] != 0)
       if(len(xmf[ind]) > 0):
           xplot = xmf[ind]
           yplot = BHSFR[i,0,ind]
           errdn = BHSFR[i,1,ind]
           errup = BHSFR[i,2,ind]
           if(i > 3):
              ax.plot(xplot,yplot[0],color=cols[i],label="z=%s" % str(c))
           else: 
               ax.plot(xplot,yplot[0],color=cols[i]) 
           if(i == 0 or i == 5 or i == 7):
              ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=cols[i], alpha=0.5, interpolate=True)
              ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=cols[i], alpha=0.5, interpolate=True)

    ax.plot([7.5,7.5], [-12.3, -7], linestyle='dotted', color='grey')
    common.prepare_legend(ax, cols[4:len(cols)], loc=3)
    plt.tight_layout()
    common.savefig(outdir, fig, 'BH-SSFR_evolution.pdf')


def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()

    read_spin = False

    if(read_spin == True):
       fields = {'galaxies': ('bh_spin', 'm_bh', 'bh_accretion_rate_hh', 'bh_accretion_rate_sb', 'sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge', 'type')}
    else:
       fields = {'galaxies': ('m_bh', 'bh_accretion_rate_hh', 'bh_accretion_rate_sb', 'sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge', 'type')}


    # Loop over redshift and subvolumes
    spinbh = np.zeros(shape = (len(zlist), 3, len(xmf)))
    spinmdot = np.zeros(shape = (len(zlist), 3, len(xmdotf))) 
    mdotmbh = np.zeros(shape = (len(zlist), 3, len(xmf)))
    mdotmbh_hh = np.zeros(shape = (len(zlist), 3, len(xmf)))
    BHSFR = np.zeros(shape = (len(zlist), 3, len(xmf)))
 
    for index, snapshot in enumerate(redshift_table[zlist]):
        print("Will read redshift", zlist[index])
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        prepare_data(hdf5_data, index, spinbh, spinmdot, mdotmbh, mdotmbh_hh, BHSFR, read_spin)

    if(read_spin == True):
       plot_spine_BH(plt, outdir, obsdir, spinbh, spinmdot, mdotmbh)
    
    plot_macc_BH(plt, outdir, obsdir, mdotmbh, mdotmbh_hh, BHSFR)


if __name__ == '__main__':
    main(*common.parse_args())
