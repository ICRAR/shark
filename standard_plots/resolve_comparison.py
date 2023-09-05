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
"""RESOLVE comparison plots"""

import functools

import numpy as np

import common
import utilities_statistics as us

##################################
# Constants
dmobs        = 0.2
HILIM        = 0.05

mlow = 5.0
mupp = 15.5
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0


# Create histogram 
def prepare_data(hdf5_data):

    (h0, volh, typeg, mdisk, mbulge, rstar_disk, m_bh, mHI, mmol_disk, mgas_disk, mHI_bulge, mmol_bulge, mgas_bulge, mhalo) = hdf5_data

    hist_bmf = np.zeros(shape = (5,len(mbins)))
    hist_bmf_sat = np.zeros(shape = (5, len(mbins)))

    hist_smf = np.zeros(shape = (5,len(mbins)))
    hist_smf_sat = np.zeros(shape = (5, len(mbins)))

    hist_himf = np.zeros(shape = (5,len(mbins)))

    Mvir_thresh = [10.0,11.4,12,13.5,16.0]
    Nbinshalo   = len(Mvir_thresh)

    mhalo  = np.log10(mhalo) - np.log10(h0)
    mbar_pseudo = mHI_bulge + mHI + mdisk + mbulge 
    #mbar_pseudo = mgas_disk + mgas_bulge + mdisk + mbulge #mHI_bulge + mHI + mdisk + mbulge
    mstar       = mdisk + mbulge

    ind = np.where((mbar_pseudo > 0) & (typeg == 0) )
    H, bins_edges = np.histogram(np.log10(mbar_pseudo[ind]) - np.log10(h0),bins=np.append(mbins,mupp))
    hist_bmf[0,:] = hist_bmf[0,:] + H

    ind = np.where((mbar_pseudo > 0) & (typeg > 0))
    H, bins_edges = np.histogram(np.log10(mbar_pseudo[ind]) - np.log10(h0),bins=np.append(mbins,mupp))
    hist_bmf_sat[0,:] = hist_bmf_sat[0,:] + H

    ind = np.where(mstar > 0)
    H, bins_edges = np.histogram(np.log10(mstar[ind]) - np.log10(h0),bins=np.append(mbins,mupp))
    hist_smf[0,:] = hist_smf[0,:] + H

    ind = np.where((mstar > 0) & (typeg > 0))
    H, bins_edges = np.histogram(np.log10(mstar[ind]) - np.log10(h0),bins=np.append(mbins,mupp))
    hist_smf_sat[0,:] = hist_smf_sat[0,:] + H

    ind = np.where(mHI_bulge+mHI > 0)
    H, bins_edges = np.histogram(np.log10(mHI_bulge[ind]+mHI[ind]) - np.log10(h0),bins=np.append(mbins,mupp))
    hist_himf[0,:] = hist_himf[0,:] + H

    sigma = 0.01
    for i in range(1,Nbinshalo):
        ind = np.where((mbar_pseudo > 0) & (mhalo > Mvir_thresh[i-1]) & (mhalo < Mvir_thresh[i]) & (typeg == 0))
        ran_err = np.random.normal(0.0, sigma, len(mbar_pseudo[ind]))
        H, bins_edges = np.histogram(np.log10(mbar_pseudo[ind]) - np.log10(h0) + ran_err,bins=np.append(mbins,mupp))
        hist_bmf[i,:] = hist_bmf[i,:] + H

        ind = np.where((mbar_pseudo > 0) & (mhalo > Mvir_thresh[i-1]) & (mhalo < Mvir_thresh[i]) & (typeg > 0))
        ran_err = np.random.normal(0.0, sigma, len(mbar_pseudo[ind]))
        H, bins_edges = np.histogram(np.log10(mbar_pseudo[ind]) - np.log10(h0) + ran_err,bins=np.append(mbins,mupp))
        hist_bmf_sat[i,:] = hist_bmf_sat[i,:] + H

        ind = np.where((mstar > 0) & (mhalo > Mvir_thresh[i-1]) & (mhalo < Mvir_thresh[i]))
        ran_err = np.random.normal(0.0, sigma, len(mbar_pseudo[ind]))
        H, bins_edges = np.histogram(np.log10(mstar[ind]) - np.log10(h0) + ran_err,bins=np.append(mbins,mupp))
        hist_smf[i,:] = hist_smf[i,:] + H

        ind = np.where((mstar > 0) & (typeg > 0) & (mhalo > Mvir_thresh[i-1]) & (mhalo < Mvir_thresh[i]))
        ran_err = np.random.normal(0.0, sigma, len(mbar_pseudo[ind]))
        H, bins_edges = np.histogram(np.log10(mstar[ind]) - np.log10(h0) + ran_err,bins=np.append(mbins,mupp))
        hist_smf_sat[i,:] = hist_smf_sat[i,:] + H

        ind = np.where((mHI_bulge+mHI > 0) & (mhalo > Mvir_thresh[i-1]) & (mhalo < Mvir_thresh[i]))
        ran_err = np.random.normal(0.0, sigma, len(mbar_pseudo[ind]))
        H, bins_edges = np.histogram(np.log10(mHI_bulge[ind]+mHI[ind]) - np.log10(h0) + ran_err,bins=np.append(mbins,mupp))
        hist_himf[i,:] = hist_himf[i,:] + H

    mHImhalo      = np.zeros(shape = (3,3,len(xmf)))
    mHIms         = np.zeros(shape = (3,3,len(xmf)))
    mHImhalo_true = np.zeros(shape = (3,3,len(xmf)))
    mHIms_true    = np.zeros(shape = (3,3,len(xmf)))

    mass    =  np.zeros(shape = (len(typeg)))
    gasfrac =  np.zeros(shape = (len(typeg)))

    ind = np.where(mdisk+mbulge > 0)
    mass[ind] = np.log10(mdisk[ind]+mbulge[ind]/h0)

    ind = np.where((mHI_bulge+mHI > 0) & (mdisk+mbulge > 0))
    gasfrac[ind] = (mHI_bulge[ind]+mHI[ind])/(mdisk[ind]+mbulge[ind])

    #apply gas fraction limit of RESOLVE 0.05:
    ind = np.where(gasfrac < HILIM)
    gasfrac[ind] = HILIM
    HIcorr = gasfrac * (mdisk+mbulge) / h0
    HItrue = (mHI_bulge+mHI) / h0

    bin_it = functools.partial(us.wmedians, xbins=xmf)

    ind = np.where((gasfrac > 0) & (mass > 8))
    mHIms[0,:] = bin_it(x=mass[ind], y=np.log10(HIcorr[ind]))
    mHImhalo[0,:] = bin_it(x=mhalo[ind], y=np.log10(HIcorr[ind]))
    mHIms_true[0,:] = bin_it(x=mass[ind], y=np.log10(HItrue[ind]))
    mHImhalo_true[0,:] = bin_it(x=mhalo[ind], y=np.log10(HItrue[ind]))

    ind = np.where((gasfrac > 0) & (typeg == 0) & (mass > 8))
    mHIms[1,:] = bin_it(x=mass[ind], y=np.log10(HIcorr[ind]))
    mHImhalo[1,:] = bin_it(x=mhalo[ind], y=np.log10(HIcorr[ind]))
    mHIms_true[1,:] = bin_it(x=mass[ind], y=np.log10(HItrue[ind]))
    mHImhalo_true[1,:] = bin_it(x=mhalo[ind], y=np.log10(HItrue[ind]))

    ind = np.where((gasfrac > 0) & (typeg > 0) & (mass >8))
    mHIms[2,:] = bin_it(x=mass[ind], y=np.log10(HIcorr[ind]))
    mHImhalo[2,:] = bin_it(x=mhalo[ind], y=np.log10(HIcorr[ind]))
    mHIms_true[2,:] = bin_it(x=mass[ind], y=np.log10(HItrue[ind]))
    mHImhalo_true[2,:] = bin_it(x=mhalo[ind], y=np.log10(HItrue[ind]))

    ETGsmhalo      = np.zeros(shape = (3,len(xmf)))
    LTGsmhalo      = np.zeros(shape = (3,len(xmf)))

    ind = np.where(mass > 8.5)
    result = us.fractions(x=mhalo[ind],y=mbulge[ind]/(mdisk[ind]+mbulge[ind]), xbins=xmf, ythresh=0.5)
    ETGsmhalo[0,:] = result
    LTGsmhalo[0,:] = 1.0-result

    ind = np.where((mass > 8.5) & (mbar_pseudo > 1e10))
    result = us.fractions(x=mhalo[ind],y=mbulge[ind]/(mdisk[ind]+mbulge[ind]), xbins=xmf, ythresh=0.5)
    ETGsmhalo[1,:] = result
    LTGsmhalo[1,:] = 1.0-result

    ind = np.where((mass > 8.5) & (mbar_pseudo < 1e10))
    result = us.fractions(x=mhalo[ind],y=mbulge[ind]/(mdisk[ind]+mbulge[ind]), xbins=xmf, ythresh=0.5)
    ETGsmhalo[2,:] = result
    LTGsmhalo[2,:] = 1.0-result

    if(volh > 0.):
        vol = volh/pow(h0,3.)  # In Mpc^3
        for i in range(0,Nbinshalo):
            hist_smf[i,:]      = hist_smf[i,:]/vol/dm
            hist_smf_sat[i,:]  = hist_smf_sat[i,:]/vol/dm
            hist_bmf[i,:]      = hist_bmf[i,:]/vol/dm
            hist_bmf_sat[i,:]  = hist_bmf_sat[i,:]/vol/dm
            hist_himf[i,:]     = hist_himf[i,:]/vol/dm

            #take logs
            ind = np.where(hist_bmf[i,:] > 0.)
            hist_bmf[i,ind] = np.log10(hist_bmf[i,ind])
            ind = np.where(hist_bmf_sat[i,:] > 0.)
            hist_bmf_sat[i,ind] = np.log10(hist_bmf_sat[i,ind])
            ind = np.where(hist_smf[i,:] > 0.)
            hist_smf[i,ind] = np.log10(hist_smf[i,ind])
            ind = np.where(hist_smf_sat[i,:] > 0.)
            hist_smf_sat[i,ind] = np.log10(hist_smf_sat[i,ind])
            ind = np.where(hist_himf[i,:] > 0.)
            hist_himf[i,ind] = np.log10(hist_himf[i,ind])

    return (hist_smf, hist_smf_sat, hist_bmf, hist_bmf_sat, mHIms, mHIms_true,
     mHImhalo, mHImhalo_true, ETGsmhalo, LTGsmhalo, hist_himf)

def _mf_obs_as_errorbar(ax, scale_factor, x, y, yerrdn, yerrup, color, marker,
                        yerrdn_val=None, **errorbar_kwargs):

    # The actual numbers that will go to the error bars
    y_plot = np.zeros(shape = (len(x)))
    yerrdn_plot = np.zeros(shape = (len(x)))
    yerrup_plot = np.zeros(shape = (len(x)))

    # Cleanup the low-end errors if necessary
    if yerrdn_val is not None:
        ind = np.where((y > 0) & (yerrdn == 0))
        if yerrdn_val is y:
            yerrdn_val = 0.01 #y[ind]
        yerrdn[ind] =  yerrdn_val

    # log and scale the rest of the values
    ind = np.where(y > 0)
    y_plot[ind] = np.log10(y[ind] / scale_factor / dmobs)
    yerrdn_plot[ind] = np.log10(yerrdn[ind] / scale_factor / dmobs)
    yerrup_plot[ind] = np.log10(yerrup[ind] / scale_factor / dmobs)

    common.errorbars(ax, x, y_plot, yerrdn_plot, yerrup_plot, color, marker,
                     condition=(y_plot != 0), **errorbar_kwargs)


def resolve_mf_obs_as_errorbar(ax, x, y, yerrdn, yerrup, color, marker,
                               yerrdn_val=None, **errorbar_kwargs):
    VOLRES = 14011.0  #Mpc^3
    _mf_obs_as_errorbar(ax, VOLRES, x, y, yerrdn, yerrup, color, marker,
                        yerrdn_val=yerrdn_val, **errorbar_kwargs)

def eco_mf_obs_as_errorbar(ax, x, y, yerrdn, yerrup, color, marker,
                           yerrdn_val=None, **errorbar_kwargs):
    VOLECO = 457956.0 #Mpc^3
    _mf_obs_as_errorbar(ax, VOLECO, x, y, yerrdn, yerrup, color, marker,
                        yerrdn_val=yerrdn_val, **errorbar_kwargs)



def _load_resolve_mf_obs(obsdir, fname, cols):
    return common.load_observation(obsdir, 'RESOLVE/massfuncs/' + fname, cols)

def _resolve_obs_as_errorbars(obsdir, ax, fname, cols, color, marker,
                              err_absolute=True, **errorbar_kwargs):
    x, y, yerrdn, yerrup = common.load_observation(obsdir, 'RESOLVE/' + fname, cols)
    common.errorbars(ax, x, y, yerrdn, yerrup, color, marker,
                     err_absolute=err_absolute, condition=(y > 0), **errorbar_kwargs)

def plot_smf_resolve(plt, outdir, obsdir, hist_smf, hist_smf_sat):

    load_resolve_mf_obs = functools.partial(_load_resolve_mf_obs, obsdir)

#   Plots global mass densities
    fig = plt.figure(figsize=(9.5,10.5))
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(\\rm dn/dM / Mpc^{-3} dex^{-1})$"
    xmin, xmax, ymin, ymax = 8, 12, -6, -1

    # all halos ##################################
    ax = fig.add_subplot(321)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=ytit)
    xleg= xmax - 0.3*(xmax-xmin)
    yleg= ymax - 0.1*(ymax-ymin)
    ax.text(xleg,yleg, '$\\rm all\, halos$')

    #Predicted SMHM
    ind = np.where(hist_smf[0,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf[0,ind]
    ax.errorbar(xplot,yplot[0],color='k')

    ind = np.where(hist_smf_sat[0,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf_sat[0,ind]
    ax.errorbar(xplot,yplot[0],color='k', linestyle="dashed")

    # RESOLVE observations
    M, No, Nodn, Noup = load_resolve_mf_obs('smassfunction_resolve.txt', [0,1,2,3])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'grey', 'o', yerrdn_val=No)

    M, No, Nodn, Noup = load_resolve_mf_obs('smassfunction_eco.txt', [0,1,2,3])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'grey', 's', yerrdn_val=No)

    # low mass halos ##################################
    ax = fig.add_subplot(322)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=None)
    xleg= xmax - 0.63*(xmax-xmin)
    yleg= ymax - 0.1*(ymax-ymin)
    ax.text(xleg,yleg, '$11<\\rm log_{10}(M_{\\rm halo}/M_{\odot})<11.4$')

    #Predicted SMHM
    ind = np.where(hist_smf[1,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf[1,ind]
    ax.errorbar(xplot,yplot[0],color='b')

    ind = np.where(hist_smf_sat[1,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf_sat[1,ind]
    ax.errorbar(xplot,yplot[0],color='b', linestyle="dashed")

    #RESOLVE observations
    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('smassfunctionlowmasshalos_resolve.txt', [0,1,2,3,7,8,9])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'b', 'o', yerrdn_val=0.1)
    resolve_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'b', 'o', yerrdn_val=0.1, fillstyle='full', markersize=3)

    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('smassfunctionlowmasshalos_eco.txt', [0,1,2,3,7,8,9])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'b', 's', yerrdn_val=0.1)
    eco_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'b', 'o', yerrdn_val=0.1, fillstyle='full', markersize=3)


    # medium mass halos ##################################
    ax = fig.add_subplot(323)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=ytit)
    ax.text(xleg, yleg, '$11.4<\\rm log_{10}(M_{\\rm halo}/M_{\odot})<12$')

    #Predicted SMHM
    ind = np.where(hist_smf[2,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf[2,ind]
    ax.errorbar(xplot,yplot[0],color='g')

    ind = np.where(hist_smf_sat[2,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf_sat[2,ind]
    ax.errorbar(xplot,yplot[0],color='g', linestyle="dashed")

    #RESOLVE observations
    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('smassfunctionmedmasshalos_resolve.txt', [0,1,2,3,7,8,9])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'g', 'o', yerrdn_val=0.1)
    resolve_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'g', 'o', fillstyle='full', markersize=3)

    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('smassfunctionmedmasshalos_eco.txt', [0,1,2,3,7,8,9])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'g', 's', yerrdn_val=0.1)
    eco_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'g', 's', fillstyle='full', markersize=3)


    # medium high mass halos ##################################
    ax = fig.add_subplot(324)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit=None)
    ax.text(xleg, yleg, '$12<\\rm log_{10}(M_{\\rm halo}/M_{\odot})<13.5$')

    #Predicted SMHM
    ind = np.where(hist_smf[3,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf[3,ind]
    ax.errorbar(xplot,yplot[0],color='r', label="all galaxies")

    ind = np.where(hist_smf_sat[3,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf_sat[3,ind]
    ax.errorbar(xplot,yplot[0],color='r', linestyle="dashed",label="satellites")

    #RESOLVE observations
    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('smassfunctionhighmasshalos_resolve.txt', [0,1,2,3,7,8,9])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'r', 'o', yerrdn_val=0.1, label="RESOLVE all")
    resolve_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'r', 'o', fillstyle='full', markersize=3, label="RESOLVE satellites")

    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('smassfunctionhighmasshalos_eco.txt', [0,1,2,3,7,8,9])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'r', 's', yerrdn_val=No, label="ECO all")
    eco_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'r', 's', fillstyle='full', markersize=3, label="ECO satellites")

    common.prepare_legend(ax, ['k','k','k','k','k','k'], loc=2, bbox_to_anchor=(0.1, -0.4))

    # high mass halos ##################################
    ax = fig.add_subplot(325)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    ax.text(xleg,yleg, '$\\rm log_{10}(M_{\\rm halo}/M_{\odot})>13.5$')

    #Predicted SMHM
    ind = np.where(hist_smf[4,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf[4,ind]
    ax.errorbar(xplot,yplot[0],color='orange')

    ind = np.where(hist_smf_sat[4,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf_sat[4,ind]
    ax.errorbar(xplot,yplot[0],color='orange', linestyle="dashed")

    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('smassfunctionclusterhalos_eco.txt', [0,1,2,3,7,8,9])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'orange', 's', yerrdn_val=No)
    eco_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'orange', 's', fillstyle='full', markersize=3)

    common.savefig(outdir, fig, "smf_resolve.pdf")


def plot_bmf_resolve(plt, outdir, obsdir, hist_bmf, hist_bmf_sat):

    load_resolve_mf_obs = functools.partial(_load_resolve_mf_obs, obsdir)

    fig = plt.figure(figsize=(9.5,10.5))

    xtit="$\\rm log_{10} (\\rm M^{\\prime}_{\\rm bar}/M_{\odot})$"
    ytit="$\\rm log_{10}(\\rm dn/dM / Mpc^{-3} dex^{-1})$"

    xmin, xmax, ymin, ymax = 9, 12, -6, -1
    xleg= xmax - 0.3 * (xmax - xmin)
    yleg= ymax - 0.1 * (ymax - ymin)

    # all halos ##################################
    ax = fig.add_subplot(321)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=ytit)
    ax.tick_params(labelsize=13)
    ax.text(xleg, yleg, '$\\rm all\, halos$')

    #Predicted SMHM
    ind = np.where(hist_bmf[0,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[0,ind]
    ax.errorbar(xplot,yplot[0],color='k')

    ind = np.where(hist_bmf_sat[0,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf_sat[0,ind]
    ax.errorbar(xplot,yplot[0],color='k', linestyle="dashed")

    #RESOLVE observations
    M, No, Nodn, Noup = load_resolve_mf_obs('bmassfunction_resolve.txt', [0,1,2,3])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'grey', 'o', yerrdn_val=No)

    M, No, Nodn, Noup = load_resolve_mf_obs('bmassfunction_eco.txt', [0,1,2,3])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'grey', 's', yerrdn_val=No)


    # low mass halos ##################################
    ax = fig.add_subplot(322)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=None)
    xleg= xmax - 0.63 * (xmax - xmin)
    yleg= ymax - 0.1 * (ymax - ymin)
    ax.text(xleg,yleg, '$11<\\rm log_{10}(M_{\\rm halo}/M_{\odot})<11.4$')

    #Predicted SMHM
    ind = np.where(hist_bmf[1,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[1,ind]
    ax.errorbar(xplot,yplot[0],color='b')
    print("Baryon MF mass bin 1")
    for a,b,c in zip(xmf, hist_bmf[1,:], hist_bmf_sat[1,:]):
        print(a,b,c)
    ind = np.where(hist_bmf_sat[1,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf_sat[1,ind]
    ax.errorbar(xplot,yplot[0],color='b', linestyle="dashed")

    #RESOLVE observations
    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionlowmasshalos_resolve.txt', [0,4,5,6,7,8,9])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'b', 'o', yerrdn_val=0.1)
    resolve_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'b', 'o', yerrdn_val=0.1, fillstyle='full', markersize=3)

    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionlowmasshalos_eco.txt', [0,1,2,3,7,8,9])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'b', 's', yerrdn_val=0.1)
    eco_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'b', 's', fillstyle='full', markersize=3)


    # medium mass halos ##################################
    ax = fig.add_subplot(323)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=ytit)
    ax.text(xleg, yleg, '$11.4<\\rm log_{10}(M_{\\rm halo}/M_{\odot})<12$')

    #Predicted SMHM
    ind = np.where(hist_bmf[2,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[2,ind]
    ax.errorbar(xplot,yplot[0],color='g')

    ind = np.where(hist_bmf_sat[2,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf_sat[2,ind]
    ax.errorbar(xplot,yplot[0],color='g', linestyle="dashed")

    print("Baryon MF mass bin 2")
    for a,b,c in zip(xmf, hist_bmf[2,:], hist_bmf_sat[2,:]):
        print(a,b,c)

    #RESOLVE observations
    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionmedmasshalos_resolve.txt', [0,1,2,3,7,8,9])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'g', 'o', yerrdn_val=0.1)
    resolve_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'g', 'o', fillstyle='full', markersize=3)

    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionmedmasshalos_eco.txt', [0,1,2,3,7,8,9])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'g', 's', yerrdn_val=0.1)
    eco_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'g', 's', fillstyle='full', markersize=3)


    # medium high mass halos ##################################
    ax = fig.add_subplot(324)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=xtit, ytit=None)
    ax.text(xleg, yleg, '$12<\\rm log_{10}(M_{\\rm halo}/M_{\odot})<13.5$')

    #Predicted SMHM
    ind = np.where(hist_bmf[3,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[3,ind]
    ax.errorbar(xplot,yplot[0],color='r', label="all galaxies")

    ind = np.where(hist_bmf_sat[3,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf_sat[3,ind]
    ax.errorbar(xplot,yplot[0],color='r', linestyle="dashed",label="satellites")
    print("Baryon MF mass bin 3")
    for a,b,c in zip(xmf, hist_bmf[3,:], hist_bmf_sat[3,:]):
        print(a,b,c)

    #RESOLVE observations
    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionhighmasshalos_resolve.txt', [0,1,2,3,7,8,9])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'r', 'o', yerrdn_val=0.1, label="RESOLVE all")
    resolve_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'r', 'o', fillstyle='full', markersize=3, label="RESOLVE satellites")

    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionhighmasshalos_eco.txt', [0,1,2,3,7,8,9])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'r', 's', yerrdn_val=No, label="ECO all")
    eco_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'r', 's', fillstyle='full', markersize=3, label="ECO satellites")

    common.prepare_legend(ax, ['k','k','k','k','k','k'], loc=2, bbox_to_anchor=(0.1, -0.4))

    # high mass halos ##################################
    ax = fig.add_subplot(325)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    ax.text(xleg, yleg, '$\\rm log_{10}(M_{\\rm halo}/M_{\odot})>13.5$')

    #Predicted SMHM
    ind = np.where(hist_bmf[4,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[4,ind]
    ax.errorbar(xplot,yplot[0],color='orange')

    ind = np.where(hist_bmf_sat[4,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf_sat[4,ind]
    ax.errorbar(xplot,yplot[0],color='orange', linestyle="dashed")

    print("Baryon MF mass bin 4")
    for a,b,c in zip(xmf, hist_bmf[4,:], hist_bmf_sat[4,:]):
        print(a,b,c)


    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionclusterhalos_eco.txt', [0,1,2,3,7,8,9])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'orange', 's', yerrdn_val=No)
    eco_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'orange', 's', fillstyle='full', markersize=3, label="ECO satellites")

    common.savefig(outdir, fig, 'bmf_resolve.pdf')


    #plot only bins in halo mass
    fig = plt.figure(figsize=(8.3,7.5))

    xtit="$\\rm log_{10}(\\rm M^{\\prime}_{\\rm bar}/M_{\odot})$"
    ytit="$\\rm log_{10}(\\phi/{\\rm dlog_{10}M^{\\prime}_{\\rm bar}}/Mpc^{-3})$"

    xmin, xmax, ymin, ymax = 9, 12.3, -6, -0.5
    xleg= xmax - 0.3 * (xmax - xmin)
    yleg= ymax - 0.1 * (ymax - ymin)


    def plot_lagos18(ax, bin_halo=1, col='k', inc_label=True, plot_central=False):
        p_pos = bin_halo * 2 - 1
        p_pos_sat = bin_halo * 2
        x,p,ps = common.load_observation(obsdir, 'Models/SharkVariations/BMF_Lagos18.dat', [0, p_pos, p_pos_sat])

        pc = p
        ind = np.where(ps != 0)
        pc[ind]  =  np.log10(10**pc[ind] - 10**ps[ind])
       
        if(plot_central == False):
           ind = np.where(p != 0)
           xplot = x[ind]
           yplot = p[ind]
           ax.errorbar(xplot,yplot,color=col, linewidth=2, linestyle='solid', alpha=0.8, label = 'Shark v1.1 (L18) all' if inc_label else None)
        else:
            ind = np.where(pc != 0)
            xplot = x[ind]
            yplot = pc[ind]
            ax.errorbar(xplot,yplot,color=col, linewidth=2, linestyle='solid', alpha=0.8, label = 'Shark v1.1 (L18) cens' if inc_label else None)
   
        ind = np.where(ps != 0)
        xplot = x[ind]
        yplot = ps[ind]
        ax.errorbar(xplot,yplot,color=col, linewidth=2, linestyle="dashed", alpha=0.8, label = 'Shark v1.1 (L18) sats' if inc_label else None)

    # low mass halos ##################################
    ax = fig.add_subplot(221)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=ytit)
    xleg= xmax - 0.7 * (xmax - xmin)
    yleg= ymax - 0.1 * (ymax - ymin)
    ax.text(xleg,yleg, '$11<\\rm log_{10}(M_{\\rm halo}/M_{\odot})<11.4$', fontsize=12)

    #Predicted SMHM
    ind = np.where(hist_bmf[1,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[1,ind]
    ax.errorbar(xplot,yplot[0],color='r', linewidth=3.5, alpha=1)

    ind = np.where(hist_bmf_sat[1,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf_sat[1,ind]
    ax.errorbar(xplot,yplot[0],color='r', linewidth=3.5,  linestyle="dashed", alpha=1)

    plot_lagos18(ax, bin_halo=1, col='k', inc_label=False, plot_central=True)
    #RESOLVE observations
    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionlowmasshalos_resolve.txt', [0,4,5,6,7,8,9])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'MediumBlue', 'o', yerrdn_val=0.1, label="ECO cens")
    resolve_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'Gold', 'o', yerrdn_val=0.1, fillstyle='full', markersize=5, label="ECO sats")

    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionlowmasshalos_eco.txt', [0,4,5,6,7,8,9])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'MediumBlue', 's', yerrdn_val=0.1, label="RESOLVE cens")
    eco_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'Gold', 's', fillstyle='full', markersize=5, label="RESOLVE sats")

    common.prepare_legend(ax, ['k','k','k','k','k','k'], fancybox=True, framealpha=0.5, bbox_to_anchor=(0.4, 0.5))

    # medium mass halos ##################################
    ax = fig.add_subplot(222)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=None)
    ax.text(xleg, yleg, '$11.4<\\rm log_{10}(M_{\\rm halo}/M_{\odot})<12$', fontsize=12)

    #Predicted SMHM
    ind = np.where(hist_bmf[2,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[2,ind]
    ax.errorbar(xplot,yplot[0],color='r', linewidth=3.5, alpha=1)

    ind = np.where(hist_bmf_sat[2,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf_sat[2,ind]
    ax.errorbar(xplot,yplot[0],color='r', linestyle="dashed", linewidth=3.5, alpha=1)

    plot_lagos18(ax, bin_halo=2, col = 'k', inc_label=False, plot_central=True)
    #RESOLVE observations
    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionmedmasshalos_resolve.txt', [0,4,5,6,7,8,9])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'MediumBlue', 'o', yerrdn_val=0.1)
    resolve_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'Gold', 'o', fillstyle='full', markersize=5)

    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionmedmasshalos_eco.txt', [0,4,5,6,7,8,9])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'MediumBlue', 's', yerrdn_val=0.1)
    eco_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'Gold', 's', fillstyle='full', markersize=5)


    # medium high mass halos ##################################
    ax = fig.add_subplot(223)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=xtit, ytit=ytit)
    ax.text(xleg, yleg, '$12<\\rm log_{10}(M_{\\rm halo}/M_{\odot})<13.5$', fontsize=12)

    #Predicted SMHM
    ind = np.where(hist_bmf[3,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[3,ind]
    ax.errorbar(xplot,yplot[0],color='r',linewidth=3.5, alpha=1)

    ind = np.where(hist_bmf_sat[3,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf_sat[3,ind]
    ax.errorbar(xplot,yplot[0],color='r', linestyle="dashed", linewidth=3.5, alpha=1)

    plot_lagos18(ax, bin_halo=3, col='k', inc_label=False, plot_central=True)
    #RESOLVE observations
    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionhighmasshalos_eco.txt', [0,4,5,6,7,8,9])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'MediumBlue', 's', yerrdn_val=No)
    eco_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'Gold', 's', fillstyle='full', markersize=5)

    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionhighmasshalos_resolve.txt', [0,4,5,6,7,8,9])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'MediumBlue', 'o', yerrdn_val=0.1)

    resolve_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'Gold', 'o', fillstyle='full', markersize=5)

    # high mass halos ##################################
    ax = fig.add_subplot(224)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=xtit, ytit=None)
    ax.text(xleg, yleg, '$\\rm log_{10}(M_{\\rm halo}/M_{\odot})>13.5$', fontsize=12)

    #Predicted SMHM
    ind = np.where(hist_bmf[4,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[4,ind]
    ax.errorbar(xplot,yplot[0],color='red', label="Shark v2.0 cens", linewidth=3.5, alpha=1)

    ind = np.where(hist_bmf_sat[4,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf_sat[4,ind]
    ax.errorbar(xplot,yplot[0],color='red', linestyle="dashed", label="Shark v2.0 sats", linewidth=3.5, alpha=1)

    plot_lagos18(ax, bin_halo=4, col='k', inc_label=True, plot_central=True)

    M, No, Nodn, Noup, Ns, Nsdn, Nsup = load_resolve_mf_obs('bmassfunctionclusterhalos_eco.txt', [0,4,5,6,7,8,9])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'MediumBlue', 's', yerrdn_val=No)
    eco_mf_obs_as_errorbar(ax, M, Ns, Nsdn, Nsup, 'Gold', 's', fillstyle='full', markersize=5)


    common.prepare_legend(ax, ['red','red','k','k'], fancybox=True, framealpha=0.5,  bbox_to_anchor=(0.22, 0.55))

    plt.tight_layout()
    common.savefig(outdir, fig, 'bmf_resolve_massbins.pdf')



def plot_mHI_mstar_resolve(plt, outdir, obsdir, mHIms, mHIms_true):

    resolve_obs_as_errorbars = functools.partial(_resolve_obs_as_errorbars, obsdir)

    fig = plt.figure(figsize=(5,9))
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(\\rm M_{\\rm HI}/M_{\odot})$"
    
    xmin, xmax, ymin, ymax = 7.5, 12, 7, 12
    xleg = xmin + 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # centrals
    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'centrals')

    #Predicted relation
    ind = np.where(mHIms[1,0,:] != 0)
    yplot = (mHIms[1,0,ind])
    errdn = (mHIms[1,1,ind])
    errup = (mHIms[1,2,ind])
    xplot = xmf[ind]

    ax.plot(xplot,yplot[0],color='k',label="G/S $>0.05$")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='grey', interpolate=True)

    ind = np.where(mHIms_true[1,0,:] != 0)
    yplot = (mHIms_true[1,0,ind])
    xplot = xmf[ind]
    ax.errorbar(xplot,yplot[0],color='b', linestyle="dashed", label="true $M_{\\rm HI}$")

    #RESOLVE observations
    resolve_obs_as_errorbars(ax, 'himassvstellarmasscentral_resolve.txt', [0,4,2,6],
                             'orange', 'o', label='RESOLVE')
    resolve_obs_as_errorbars(ax, 'himassvstellarmasscentral_eco.txt', [0,4,2,6],
                             'orange', 's', label='ECO')

    common.prepare_legend(ax, ['k','b','orange','orange'], loc=4)

    # satellites
    ax = fig.add_subplot(212)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'satellites')

    #Predicted relation
    ind = np.where(mHIms[2,0,:] != 0)
    yplot = (mHIms[2,0,ind])
    errdn = (mHIms[2,1,ind])
    errup = (mHIms[2,2,ind])
    xplot = xmf[ind]

    ax.plot(xplot,yplot[0],color='k')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='grey', interpolate=True)

    ind = np.where(mHIms_true[2,0,:] != 0)
    yplot = (mHIms_true[2,0,ind])
    xplot = xmf[ind]
    ax.errorbar(xplot,yplot[0],color='b', linestyle="dashed")

    # RESOLVE observations
    resolve_obs_as_errorbars(ax, 'himassvstellarmasssatellite_resolve.txt', [0,4,2,6], 'orange', 'o')
    resolve_obs_as_errorbars(ax, 'himassvstellarmasssatellite_eco.txt', [0,4,2,6], 'orange', 's')

    common.savefig(outdir, fig, 'mHI-mstar_resolve.pdf')


def plot_mHI_mhalo_resolve(plt, outdir, obsdir, mHImhalo, mHImhalo_true):

    resolve_obs_as_errorbars = functools.partial(_resolve_obs_as_errorbars, obsdir)

    #   Plots gas metallicity vs. stellar mass
    fig = plt.figure(figsize=(5,9))
    xtit="$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit="$\\rm log_{10}(\\rm M_{\\rm HI}/M_{\odot})$"

    xmin, xmax, ymin, ymax = 10.5, 15.2, 7, 12
    xleg= xmin + 0.2 * (xmax - xmin)
    yleg= ymax - 0.1 * (ymax - ymin)

    # centrals
    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'centrals')

    #Predicted relation
    ind = np.where(mHImhalo[1,0,:] != 0)
    yplot = (mHImhalo[1,0,ind])
    errdn = (mHImhalo[1,1,ind])
    errup = (mHImhalo[1,2,ind])
    xplot = xmf[ind]

    ax.plot(xplot,yplot[0],color='k',label="G/S $>0.05$")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='grey', interpolate=True)

    ind = np.where(mHImhalo_true[1,0,:] != 0)
    yplot = (mHImhalo_true[1,0,ind])
    xplot = xmf[ind]
    ax.errorbar(xplot,yplot[0],color='b', linestyle="dashed", label="true $M_{\\rm HI}$")

    #RESOLVE observations
    resolve_obs_as_errorbars(ax, 'himassvgroupmasscentral_resolve.txt', [0,4,2,6], 'orange', 'o', label='RESOLVE')
    resolve_obs_as_errorbars(ax, 'himassvgroupmasscentral_eco.txt', [0,4,2,6], 'orange', 'o', label='ECO')

    common.prepare_legend(ax, ['k','b','orange','orange'], loc=4)

    # satellites
    ax = fig.add_subplot(212)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'satellites')

    #Predicted relation
    ind = np.where(mHImhalo[2,0,:] != 0)
    yplot = (mHImhalo[2,0,ind])
    errdn = (mHImhalo[2,1,ind])
    errup = (mHImhalo[2,2,ind])
    xplot = xmf[ind]

    ax.plot(xplot,yplot[0],color='k')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='grey', interpolate=True)

    ind = np.where(mHImhalo_true[2,0,:] != 0)
    yplot = (mHImhalo_true[2,0,ind])
    xplot = xmf[ind]
    ax.errorbar(xplot,yplot[0],color='b', linestyle="dashed")

    #RESOLVE observations
    resolve_obs_as_errorbars(ax, 'himassvgroupmasssatellite_resolve.txt', [0,4,2,6], 'orange', 'o')
    resolve_obs_as_errorbars(ax, 'himassvgroupmasssatellite_eco.txt', [0,4,2,6], 'orange', 's')

    #Save figure
    common.savefig(outdir, fig, 'mHI-mhalo_resolve.pdf')


def plot_bt_resolve(plt, outdir, obsdir, ETGsmhalo, LTGsmhalo):

    resolve_obs_as_errorbars = functools.partial(_resolve_obs_as_errorbars, obsdir)

    fig = plt.figure(figsize=(5,9))
    xtit="$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit="$\\rm frequency$"
    xmin, xmax, ymin, ymax = 11.5, 15, -0.05, 1.05
    xleg= xmax - 0.5 * (xmax - xmin)
    yleg= ymax - 0.1 * (ymax - ymin)

    # LTG ##################################
    ax = fig.add_subplot(311)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, '$\\rm all\, masses$')

    #Predicted size-mass for disks
    ind = np.where(ETGsmhalo[0,:] >= 0)
    xplot = xmf[ind]
    yplot = ETGsmhalo[0,ind]
    ax.plot(xplot,yplot[0],'r', label ='SHArk ETGs')

    ind = np.where(LTGsmhalo[0,:] >= 0)
    xplot = xmf[ind]
    yplot = LTGsmhalo[0,ind]
    ax.plot(xplot,yplot[0],'b', linestyle='dashed',label ='SHArk LTGs')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    resolve_obs_as_errorbars(ax, 'ETall_frac.txt', [0, 1, 2, 3], 'r', 'o', err_absolute=False, label='RESOLVE/ECO ETGs')
    resolve_obs_as_errorbars(ax, 'LTall_frac.txt', [0, 1, 2, 3], 'b', 's', err_absolute=False, label='RESOLVE/ECO LTGs')

    common.prepare_legend(ax, ['r','b','r','b'], loc=2, bbox_to_anchor=(0.0, 1.45))

    # LTG ##################################
    ax = fig.add_subplot(312)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, '$M^{\\prime}_{\\rm bar}/M_{\odot} > 10^{10}$')

    #Predicted size-mass for disks
    ind = np.where(ETGsmhalo[1,:] >= 0)
    xplot = xmf[ind]
    yplot = ETGsmhalo[1,ind]
    ax.plot(xplot,yplot[0],'r')

    ind = np.where(LTGsmhalo[1,:] >= 0)
    xplot = xmf[ind]
    yplot = LTGsmhalo[1,ind]
    ax.plot(xplot,yplot[0],'b', linestyle='dashed')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    resolve_obs_as_errorbars(ax, 'EThimbary_frac.txt', [0, 1, 2, 3], 'r', 'o', err_absolute=False)
    resolve_obs_as_errorbars(ax, 'LThimbary_frac.txt', [0, 1, 2, 3], 'b', 's', err_absolute=False)

    # LTG ##################################
    ax = fig.add_subplot(313)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, '$M^{\\prime}_{\\rm bar}/M_{\odot} < 10^{10}$')

    #Predicted size-mass for disks
    ind = np.where(ETGsmhalo[2,:] >= 0)
    xplot = xmf[ind]
    yplot = ETGsmhalo[2,ind]
    ax.plot(xplot,yplot[0],'r')

    ind = np.where(LTGsmhalo[2,:] >= 0)
    xplot = xmf[ind]
    yplot = LTGsmhalo[2,ind]
    ax.plot(xplot,yplot[0],'b', linestyle='dashed')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    resolve_obs_as_errorbars(ax, 'ETlowmbary_frac.txt', [0, 1, 2, 3], 'r', 'o', err_absolute=False)
    resolve_obs_as_errorbars(ax, 'LTlowmbary_frac.txt', [0, 1, 2, 3], 'b', 's', err_absolute=False)

    common.savefig(outdir, fig, 'bt_resolve.pdf')


def plot_bmf_resolve_baryon_components(plt, outdir, obsdir, hist_bmf, hist_smf, hist_himf):

    load_resolve_mf_obs = functools.partial(_load_resolve_mf_obs, obsdir)

    fig = plt.figure(figsize=(9.5,10.5))
    xtit = "$\\rm log_{10} (\\rm M^{\\prime}_{\\rm bar}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\\rm dn/dM / Mpc^{-3} dex^{-1})$"
    xmin, xmax, ymin, ymax = 9, 12, -6, -1
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # all halos ##################################
    ax = fig.add_subplot(321)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=ytit)
    ax.text(xleg, yleg, '$\\rm all\, halos$')

    #Predicted SMHM
    ind = np.where(hist_bmf[0,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[0,ind]
    ax.errorbar(xplot,yplot[0],color='k')

    ind = np.where(hist_smf[0,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf[0,ind]
    ax.errorbar(xplot,yplot[0],color='k', linestyle="dotted")

    ind = np.where(hist_himf[0,:] != 0)
    xplot = xmf[ind]
    yplot = hist_himf[0,ind]
    ax.errorbar(xplot,yplot[0],color='k', linestyle="dashed")

    # RESOLVE observations
    M, No, Nodn, Noup = load_resolve_mf_obs('bmassfunction_resolve.txt', [0,1,2,3])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'grey', 'o', yerrdn_val=No)
    M, No, Nodn, Noup = load_resolve_mf_obs('bmassfunction_eco.txt', [0,1,2,3])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'grey', 's', yerrdn_val=No)

    # low mass halos ##################################
    ax = fig.add_subplot(322)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=None)
    xleg = xmax - 0.63 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)
    ax.text(xleg, yleg, '$11<\\rm log_{10}(M_{\\rm halo}/M_{\odot})<11.4$')

    #Predicted SMHM
    ind = np.where(hist_bmf[1,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[1,ind]
    ax.errorbar(xplot,yplot[0],color='b')

    ind = np.where(hist_smf[1,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf[1,ind]
    ax.errorbar(xplot,yplot[0],color='b', linestyle="dotted")

    ind = np.where(hist_himf[1,:] != 0)
    xplot = xmf[ind]
    yplot = hist_himf[1,ind]
    ax.errorbar(xplot,yplot[0],color='b', linestyle="dashed")

    #RESOLVE observations
    M, No, Nodn, Noup = load_resolve_mf_obs('bmassfunctionlowmasshalos_resolve.txt', [0,1,2,3])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'b', 'o', yerrdn_val=0.1)
    M, No, Nodn, Noup = load_resolve_mf_obs('bmassfunctionlowmasshalos_eco.txt', [0,1,2,3])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'b', 's', yerrdn_val=0.1)

    # medium mass halos ##################################
    ax = fig.add_subplot(323)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit=None, ytit=ytit)
    ax.tick_params(labelsize=13)
    ax.text(xleg, yleg, '$11.4<\\rm log_{10}(M_{\\rm halo}/M_{\odot})<12$')

    #Predicted SMHM
    ind = np.where(hist_bmf[2,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[2,ind]
    ax.errorbar(xplot,yplot[0],color='g')

    ind = np.where(hist_smf[2,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf[2,ind]
    ax.errorbar(xplot,yplot[0],color='g', linestyle="dotted")

    ind = np.where(hist_himf[2,:] != 0)
    xplot = xmf[ind]
    yplot = hist_himf[2,ind]
    ax.errorbar(xplot,yplot[0],color='g', linestyle="dashed")

    #RESOLVE observations
    M, No, Nodn, Noup = load_resolve_mf_obs('bmassfunctionmedmasshalos_resolve.txt', [0,1,2,3])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'g', 'o', yerrdn_val=0.1)
    M, No, Nodn, Noup = load_resolve_mf_obs('bmassfunctionmedmasshalos_eco.txt', [0,1,2,3])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'g', 's', yerrdn_val=0.1)


    # medium high mass halos ##################################
    ax = fig.add_subplot(324)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit=None)
    ax.text(xleg, yleg, '$12<\\rm log_{10}(M_{\\rm halo}/M_{\odot})<13.5$')

    #Predicted SMHM
    ind = np.where(hist_bmf[3,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[3,ind]
    ax.errorbar(xplot,yplot[0],color='r', label="all galaxies")

    ind = np.where(hist_smf[3,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf[3,ind]
    ax.errorbar(xplot,yplot[0],color='r', linestyle="dotted",label="stellar mass")

    ind = np.where(hist_himf[3,:] != 0)
    xplot = xmf[ind]
    yplot = hist_himf[3,ind]
    ax.errorbar(xplot,yplot[0],color='r', linestyle="dashed",label='HI mass')

    #RESOLVE observations
    M, No, Nodn, Noup = load_resolve_mf_obs('bmassfunctionhighmasshalos_resolve.txt', [0,1,2,3])
    resolve_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'r', 'o', yerrdn_val=0.1, label="RESOLVE all")
    M, No, Nodn, Noup = load_resolve_mf_obs('bmassfunctionhighmasshalos_eco.txt', [0,1,2,3])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'r', 's', yerrdn_val=No, label="ECO all")

    common.prepare_legend(ax, ['k','k','k','k','k','k'], loc=2, bbox_to_anchor=(0.1, -0.4))


    # high mass halos ##################################
    ax = fig.add_subplot(325)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    ax.text(xleg, yleg, '$\\rm log_{10}(M_{\\rm halo}/M_{\odot})>13.5$')

    #Predicted SMHM
    ind = np.where(hist_bmf[4,:] != 0)
    xplot = xmf[ind]
    yplot = hist_bmf[4,ind]
    ax.errorbar(xplot,yplot[0],color='orange')

    ind = np.where(hist_smf[4,:] != 0)
    xplot = xmf[ind]
    yplot = hist_smf[4,ind]
    ax.errorbar(xplot,yplot[0],color='orange', linestyle="dotted")

    ind = np.where(hist_himf[4,:] != 0)
    xplot = xmf[ind]
    yplot = hist_himf[4,ind]
    ax.errorbar(xplot,yplot[0],color='orange', linestyle="dashed")

    M, No, Nodn, Noup = load_resolve_mf_obs('bmassfunctionclusterhalos_eco.txt', [0,1,2,3])
    eco_mf_obs_as_errorbar(ax, M, No, Nodn, Noup, 'orange', 's', yerrdn_val=No, label="ECO all")

    common.savefig(outdir, fig, 'bmf_resolve_baryon-components.pdf')


def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge', 'rstar_disk',
                           'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge',
                           'mvir_hosthalo')}

    hdf5_data = common.read_data(modeldir, redshift_table[0], fields, subvols)

    (hist_smf, hist_smf_sat, hist_bmf, hist_bmf_sat,
     mHIms, mHIms_true, mHImhalo, mHImhalo_true,
     ETGsmhalo, LTGsmhalo, hist_himf) = prepare_data(hdf5_data)

    plot_smf_resolve(plt, outdir, obsdir, hist_smf, hist_smf_sat)
    plot_bmf_resolve(plt, outdir, obsdir, hist_bmf, hist_bmf_sat)
    plot_mHI_mstar_resolve(plt, outdir, obsdir, mHIms, mHIms_true)
    plot_mHI_mhalo_resolve(plt, outdir, obsdir, mHImhalo, mHImhalo_true)
    plot_bt_resolve(plt, outdir, obsdir, ETGsmhalo, LTGsmhalo)
    plot_bmf_resolve_baryon_components(plt, outdir, obsdir, hist_bmf, hist_smf, hist_himf)

if __name__ == '__main__':
    main(*common.parse_args())
