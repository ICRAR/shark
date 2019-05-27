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

import common
import utilities_statistics as us

##################################
# Constants
mlow = 8.0
mupp = 12.0
dm = 1.0
mbins = np.arange(mlow, mupp, dm)
xmf = mbins + dm/2.0

rlow = 0.0
rupp = 5.0
dr = 0.5
rbins = np.arange(rlow, rupp, dr)
xrf = rbins + dr/2.0

GyrtoYr  = 1e9
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun

def add_observations_to_plot(obsdir, fname, ax, marker, label, color='k', err_absolute=False):
    fname = '%s/Gas/%s' % (obsdir, fname)
    x, y, yerr_down, yerr_up = common.load_observation(obsdir, fname, (0, 1, 2, 3))
    common.errorbars(ax, x, y, yerr_down, yerr_up, color, marker, label=label, err_absolute=err_absolute)

def prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit):
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    xleg = xmax - 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)
    #ax.text(xleg, yleg, 'z=0')

def prepare_data(hdf5_data, fradii, index):

    bin_it = functools.partial(us.wmedians, xbins=xmf)
    stack  = functools.partial(us.stacking, xbins=xmf)

    # Unpack data
    (h0, _, typeg, mdisk, mbulge, _, _, mHI, mH2, mgas,
     mHI_bulge, mH2_bulge, mgas_bulge, mvir, sfrd, sfrb, 
     x, y, z, vvir) = hdf5_data

    XH = 0.72
    h0log = np.log10(float(h0))

    rvir  = G * mvir / pow(vvir,2.0) / h0

    mstar_tot = (mdisk + mbulge) / h0
    print max(mstar_tot)
    sfr_tot = (sfrd + sfrb) / h0 / GyrtoYr

    #define main sequence first
    inms = np.where((mstar_tot > 5e8) & (mstar_tot < 7e9) & (typeg == 0) & (sfr_tot > 0))
    ms = np.polyfit(np.log10(mstar_tot[inms]), np.log10(sfr_tot[inms]), 2)
    gasfracms = np.polyfit(np.log10(mstar_tot[inms]), np.log10(mgas[inms]+mgas_bulge[inms])-h0log, 2)

    print 'ms',ms
    indcen = np.where((mvir/h0 > 1e14) & (typeg == 0))
    x_cen = x[indcen]
    y_cen = y[indcen]
    z_cen = z[indcen]
    rvir_cen = rvir[indcen]

    print 'number of clusters', len(x_cen)
    nradii_this_z = np.zeros(shape = (3, len(xmf), len(xrf), len(x_cen)))
    #xy projection
    for g in range(0,len(x_cen)):
        d_all = np.sqrt(pow(x - x_cen[g], 2.0) + pow(y - y_cen[g], 2.0))/h0/rvir_cen[g]
        for i in range(0, len(xmf)):
            ind   = np.where((np.log10(mstar_tot) >= xmf[i] - dm/2.0) 
                    & (np.log10(mstar_tot) < xmf[i] +  dm/2.0) 
                    & (d_all < 5.5))
            #print 'number of neighbours', len(sfr_tot[ind])
            mstars_galsin   = np.log10(mstar_tot[ind])
            sfr_tot_galsin  = sfr_tot[ind]
            mgas_tot_galsin = (mgas[ind] + mgas_bulge[ind])/h0
            dist_to_ms      = sfr_tot_galsin / pow(10.0, (ms[0] * mstars_galsin**2.0 + ms[1] * mstars_galsin + ms[2]))
            dist_to_gf      = mgas_tot_galsin / pow(10.0, (gasfracms[0] * mstars_galsin**2.0 + gasfracms[1] * mstars_galsin + gasfracms[2]))
            dist_proj       = d_all[ind]
            for j in range(0, len(xrf)):
                inr = np.where((dist_proj >= xrf[j] - dr/2.0) & (dist_proj < xrf[j] + dr/2.0))
                nradii_this_z[0,i,j,g] = len(dist_proj[inr])
                inr = np.where((dist_proj >= xrf[j] - dr/2.0) & (dist_proj < xrf[j] + dr/2.0) & (dist_to_ms > 0.33))
                nradii_this_z[1,i,j,g] = len(dist_proj[inr])
                inr = np.where((dist_proj >= xrf[j] - dr/2.0) & (dist_proj < xrf[j] + dr/2.0) & (dist_to_gf > 0.33))
                nradii_this_z[2,i,j,g] = len(dist_proj[inr])

    for i in range(0, len(xmf)):
        for j in range(0, len(xrf)):
            selec_cl = np.where(nradii_this_z[0,i,j,:] > 0)
            fradii[0,0,index,i,j] = np.median(nradii_this_z[1,i,j,selec_cl] / nradii_this_z[0,i,j,selec_cl])
            fradii[0,1,index,i,j] = np.std(nradii_this_z[1,i,j,selec_cl] / nradii_this_z[0,i,j,selec_cl])
            fradii[1,0,index,i,j] = np.median(nradii_this_z[2,i,j,selec_cl] / nradii_this_z[0,i,j,selec_cl])
            fradii[1,1,index,i,j] = np.std(nradii_this_z[2,i,j,selec_cl] / nradii_this_z[0,i,j,selec_cl])

def plot_fractions_radii(plt, output_dir, fradii):

    ###################################
    #   Plots global mass densities
    fig = plt.figure(figsize=(6,7))

    plt.subplots_adjust(bottom=0.15, left=0.15)
    subplots = (321, 322, 323, 324, 325, 326)
    zs = (0, 0.3, 0.5)
    cols = ('r','g','b')
    lines = ('dotted', 'dashed', 'solid')
    xmin, xmax, ymin, ymax = 0, 5, -0.05, 1.05

    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymin + 0.1 * (ymax - ymin)

    xtitle = '$\\rm d_{\\rm proj}/cMpc$'
    ytitle = '$\\rm fraction$'

    labels = ('Main sequence', 'Gas rich')

    p = 0
    for i in range(0, len(xmf)-1):
        for j in range(0,2):
            print i,j
            ax = fig.add_subplot(subplots[p])
            if(p <= 1):
               ax.text(1,1.1,labels[j])

            if (p >= 4):
                xtit = xtitle
            else:
                xtit = ''
            if (p == 0 or p == 2 or p == 4):
                ytit = ytitle
            else:
                ytit = ''
            common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1))
            ax.text(xleg, yleg, '$M_{\\star}$=%s' % str(xmf[i]))
        
            #predicted fraction
            for z in range (0,3):
                ind = np.where(fradii[j,0,z,i,:] > 0)
                xplot = xrf[ind]
                yplot = fradii[j,0,z,i,ind]
                err   = fradii[j,1,z,i,ind]
                ax.plot(xplot, yplot[0], color=cols[z], linestyle=lines[z])
                ax.fill_between(xplot,yplot[0],yplot[0]-err[0], facecolor=cols[z], alpha=0.5,interpolate=True)
                ax.fill_between(xplot,yplot[0],yplot[0]+err[0], facecolor=cols[z], alpha=0.5,interpolate=True)

            p = p + 1


    #common.prepare_legend(ax, ['k','b','r','grey','grey'])
    common.savefig(output_dir, fig, "cluster_fractions.pdf")

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):


    plt = common.load_matplotlib()

    zlist = (0, 0.3, 0.5)

    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge', 'mvir_hosthalo', 'sfr_disk',
                           'sfr_burst', 'position_x', 'position_y', 'position_z', 'vvir_hosthalo')}

    fradii = np.zeros(shape = (2, 2, len(zlist), len(xmf), len(xrf)))
    fradii[:] = -1

    for index, snapshot in enumerate(redshift_table[zlist]):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        prepare_data(hdf5_data, fradii, index)

    plot_fractions_radii(plt, output_dir, fradii)

if __name__ == '__main__':
    main(*common.parse_args())
