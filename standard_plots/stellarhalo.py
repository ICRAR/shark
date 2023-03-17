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

import numpy as np
import collections
import logging

import common
import utilities_statistics as us


observation = collections.namedtuple('observation', 'label x y yerrup yerrdn err_absolute')

logger = logging.getLogger(__name__)

mlow = 8
mupp = 15
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

mrlow = -5
mrupp = 5
dmr = 0.2
mrbins = np.arange(mrlow,mrupp,dmr)
xmrf = mrbins + dmr/2.0

def prepare_data(hdf5_data, index, massstellarh, massstellarh_v2, hist_smf, hist_minfall, frac_sat):

    Omegab = 0.0491
    OmegaM = 0.3121

    (h0, volh, mdisk, mbulge, ms_halo, mhalo, typeg, ms_tidally_stripped, 
    id_halo, mvir_infall, mvir_subhalo, id_subhalo) = hdf5_data

    vol = volh/pow(h0,3.)  # In Mpc^3

    ids_halo = np.unique(id_halo)

    ms_gals = np.zeros(shape = (3,len(ids_halo)))
    n_gals = np.zeros(shape = (len(ids_halo)))
    #for i,g in enumerate(ids_halo):
    #    ind = np.where((id_halo == g) & (mdisk+mbulge > 0))
    #    indcen = np.where((id_halo == g) & (typeg == 0))
    #    if((np.count_nonzero(ind) > 0) & (np.count_nonzero(indcen) > 0)):
    #       ms_gals[0,i] = np.sum(mdisk[ind]+mbulge[ind])
    #       ms_gals[0,i] = ms_gals[0,i] + ms_halo[indcen]
    #       ms_gals[1,i] = ms_halo[indcen]
    #       ms_gals[2,i] = mhalo[indcen]

    i = 0
    j = 0
    for i in range(0,len(ids_halo)):
        ids = id_halo[j]
        while (id_halo[j] == ids):
            ms_gals[0,i] = ms_gals[0,i] + mdisk[j] + mbulge[j]
            n_gals[i] = n_gals[i] + 1
            if(typeg[j] == 0):
                ms_gals[1,i] = ms_halo[j]
                ms_gals[2,i] = mhalo[j]
                n_gals[i] = n_gals[i] - 1
            j = j + 1
            if(j == len(id_halo)):
               break
        i = i + 1

    for i, m in enumerate(xmf):
        ind = np.where((ms_gals[2,:] >= 10**(m - dm/2.0)) & (ms_gals[2,:] < 10**(m + dm/2.0)))
        if(len(n_gals[ind]) == 0):
            frac_sat[index,i] = -10 
        else:
            frac_sat[index,i] = np.median(n_gals[ind])
    ind = np.where(frac_sat[index,:] == 0)
    frac_sat[index,ind] = 0.01


    ind = np.where((typeg <= 0) & (ms_halo > 0))
    massstellarh[index,:] = us.wmedians(x=np.log10(mhalo[ind]) - np.log10(float(h0)),
                                   y=np.log10(ms_halo[ind]) - np.log10(mhalo[ind]),
                                   xbins=xmf)
    ind = np.where(ms_gals[1,:] == 0)
    ms_gals[1,ind] = 1.0

    ind = np.where(ms_gals[0,:] > 0)
    massstellarh_v2[index,:] = us.wmedians(x=np.log10(ms_gals[2,ind]) - np.log10(float(h0)),
                                   y=np.log10(ms_gals[1,ind]) - np.log10(ms_gals[0,ind]),
                                   xbins=xmf)

    mass = mdisk + mbulge
    ind = np.where(mass > 0)
    H, _ = np.histogram(np.log10(mass[ind]/h0),bins=np.append(mbins,mupp))
    hist_smf[index,0,:] = hist_smf[index,0,:] + H
    ind = np.where((mass > 0) & (typeg == 0))
    H, _ = np.histogram(np.log10(mass[ind]/h0),bins=np.append(mbins,mupp))
    hist_smf[index,1,:] = hist_smf[index,1,:] + H
    ind = np.where((mass > 0) & (typeg == 1))
    H, _ = np.histogram(np.log10(mass[ind]/h0),bins=np.append(mbins,mupp))
    hist_smf[index,2,:] = hist_smf[index,2,:] + H
    ind = np.where((mass > 0) & (typeg == 2))
    H, _ = np.histogram(np.log10(mass[ind]/h0),bins=np.append(mbins,mupp))
    hist_smf[index,3,:] = hist_smf[index,3,:] + H

    ind = np.where(hist_smf > 0)
    hist_smf[ind] = np.log10(hist_smf[ind]/vol/dm)

    ind = np.where(mvir_infall > 0)
    H , _ = np.histogram(np.log10(mvir_subhalo[ind]/mvir_infall[ind]), bins=np.append(mrbins,mrupp))
    hist_minfall[index,:] = hist_minfall[index,:] + H
    hist_minfall[index,:] = hist_minfall[index,:] / np.sum(hist_minfall[index,:] * dm)

    return (h0)

def load_smf_observations(obsdir, h0):

    # Wright et al. (2017, z=0). Chabrier IMF
    z0obs = []
    lm, p, dpdn, dpup = common.load_observation(obsdir, 'mf/SMF/GAMAII_BBD_GSMFs.dat', [0,1,2,3])
    hobs = 0.7
    xobs = lm + 2.0 * np.log10(hobs/h0)
    indx = np.where(p > 0)
    yobs = np.log10(p[indx]) - 3.0 * np.log10(hobs/h0)
    ydn = yobs - np.log10(p[indx]-dpdn[indx])
    yup = np.log10(p[indx]+dpup[indx]) - yobs
    z0obs.append((observation("Wright+2017", xobs[indx], yobs, ydn, yup, err_absolute=False), 'o'))

    lm, p, dpdn, dpup = common.load_observation(obsdir, 'mf/SMF/SMF_Bernardi2013_SerExp.data', [0,1,2,3])
    xobs = lm + 2.0 * np.log10(hobs/h0)
    indx = np.where(p > 0)
    yobs = np.log10(p[indx]) - 3.0 * np.log10(hobs/h0)
    ydn = yobs - np.log10(p[indx]-dpdn[indx])
    yup = np.log10(p[indx]+dpup[indx]) - yobs
    z0obs.append((observation("Bernardi+2013", xobs[indx], yobs, ydn, yup, err_absolute=False), 's'))

    return z0obs

def plot_stellarmf_z(plt, outdir, obsdir, h0, hist_smf, hist_minfall, zlist):

    (z0obs) = load_smf_observations(obsdir, h0)

    fig = plt.figure(figsize=(6,5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm M_{\\star}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 8, 13, -6, -1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = [111]
    indeces = [0]
    zs = [0]
    observations = [z0obs]

    for subplot, idx, z, obs_and_markers in zip(subplots, indeces, zs, observations):

        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(z))

        # Observations
        for obs, marker in obs_and_markers:
            common.errorbars(ax, obs.x, obs.y, obs.yerrdn, obs.yerrup, 'grey',
                             marker, err_absolute=obs.err_absolute, label=obs.label)

        # Predicted SMF
        y = hist_smf[idx,0,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'r', label='all galaxies' if idx == 0 else None)
        y = hist_smf[idx,1,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'b', linestyle='dotted', label ='type=0' if idx == 0 else None)
        y = hist_smf[idx,2,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'g', linestyle='dashed', label ='type=1' if idx == 0 else None)
        y = hist_smf[idx,3,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'DarkTurquoise', linestyle='dashdot', label ='type=2' if idx == 0 else None)

        colors = []
        if idx == 0:
            colors = ['r','b','g','DarkTurquoise']
        colors += ['grey', 'grey','grey']

        common.prepare_legend(ax, colors)

    common.savefig(outdir, fig, 'stellarmf_z_satellites_contribution.pdf')


    fig = plt.figure(figsize=(6,5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm subhalo}/M_{\\rm infall})$"
    ytit = "$\\rm PDF$"
    xmin, xmax, ymin, ymax = -5, 2, 0, 3
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = [111]
    indeces = [0]
    zs = [0]

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    cols = ['DarkRed', 'DarkOrange', 'LimeGreen', 'Turquoise', 'SteelBlue', 'Purple']
    #Predicted stellar-halo mass relation
    for i,z in enumerate(zlist):
        y = hist_minfall[i,:]
        ind = np.where(y != 0.)
        ax.plot(xmrf[ind],y[ind], linestyle='solid', color=cols[i], label = 'z=%s' % str(z))
    xp = [0,0]
    yp = [0,3]
    ax.plot(xp,yp,color='k', linestyle='dotted')

    common.prepare_legend(ax, cols)

    common.savefig(outdir, fig, 'infall_masses_distribution.pdf')

def plot_stellar_halo_z(plt, outdir, massstellarh, massstellarh_v2, zlist, frac_sat):

    fig = plt.figure(figsize=(7, 9))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo, DM}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm M_{\\star,halo}/M_{\\star,total})$"
    xmin, xmax, ymin, ymax = 10.5, 15, -4, 0

    # z=0 ##################################
    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    plt.subplots_adjust(left=0.2)

    def compute_points(massstellarh):
        ind = np.where(massstellarh[i,0,:] != 0)
        xplot = xmf[ind]
        yplot = massstellarh[i,0,ind]
        errdn = massstellarh[i,1,ind]
        errup = massstellarh[i,2,ind]
        return (xplot, yplot, errdn, errup)
    
    cols = ['DarkRed', 'DarkOrange', 'LimeGreen', 'Turquoise', 'SteelBlue', 'Purple']
    #Predicted stellar-halo mass relation
    for i,z in enumerate(zlist):
        (xplot, yplot, errdn, errup) = compute_points(massstellarh_v2)
        ax.plot(xplot, yplot[0], linestyle = 'solid', color=cols[i], label = 'z=%s' % str(z))
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=cols[i], alpha=0.5, interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=cols[i], alpha=0.5, interpolate=True)

#        (xplot, yplot, errdn, errup) = compute_points(massstellarh_v2)
#        ax.plot(xplot, yplot[0], linestyle = 'solid', color=cols[i], label = 'z=%s' % str(z))
#        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=cols[i], alpha=0.5, interpolate=True)
#        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=cols[i], alpha=0.5, interpolate=True)

    common.prepare_legend(ax, cols, loc=2)


    ax = fig.add_subplot(212)
    ytit = "$\\rm log_{10}(\\langle N_{sat}\\rangle)$"
    xmin, xmax, ymin, ymax = 10.5, 15, -0.2, 5
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    plt.subplots_adjust(left=0.2)

    for i,z in enumerate(zlist):
        print(frac_sat[i,:])
        ind = np.where(frac_sat[i,:] > -10)
        xplot = xmf[ind]
        yplot = np.log10(frac_sat[i,ind])
        ax.plot(xplot, yplot[0], linestyle = 'solid', color=cols[i], label = 'z=%s' % str(z))

    common.savefig(outdir, fig, 'stellar_halo_mass_z.pdf')


def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstellar_halo', 'mvir_hosthalo',
                           'type','mstars_tidally_stripped','id_halo_tree', 'mvir_infall_subhalo',
                           'mvir_subhalo', 'id_subhalo_tree')}

    zlist = [0, 0.5, 1] #, 2, 3, 4]
    snapshots = redshift_table[zlist]
    massstellarh    = np.zeros(shape = (len(zlist), 3, len(xmf)))
    massstellarh_v2 = np.zeros(shape = (len(zlist), 3, len(xmf)))
    hist_smf        = np.zeros(shape = (len(zlist), 4, len(mbins)))
    hist_minfall    = np.zeros(shape = (len(zlist), len(mrbins)))
    frac_sat        =  np.zeros(shape = (len(zlist), len(xmf)))

    for idx, snapshot in enumerate(snapshots):
        print("Will read and process redshift", zlist[idx])
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        (h0) = prepare_data(hdf5_data, idx, massstellarh, massstellarh_v2, hist_smf, hist_minfall, frac_sat)


    plot_stellar_halo_z(plt, outdir, massstellarh, massstellarh_v2, zlist, frac_sat)
    plot_stellarmf_z(plt, outdir, obsdir, h0, hist_smf, hist_minfall, zlist)


if __name__ == '__main__':
    main(*common.parse_args(requires_observations=True))
