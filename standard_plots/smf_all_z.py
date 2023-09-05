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

##################################
# Mass function initialization
mlow = 7
mupp = 15
dm = 0.125
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

mlow2 = 8
mupp2 = 12
dm2 = 0.15
mbins2 = np.arange(mlow2,mupp2,dm2)
xmf2 = mbins2 + dm2/2.0

def load_smf_passive_observations(obsdir, h0):

    # Weaver et al. (2022; COSMOS2020)
    def add_weaver22_data(zappend, file_read='0.2z0.5'):
        lm, pD, dn, du = np.loadtxt(obsdir+'/mf/SMF/COSMOS2020/SMF_Farmer_v2.1_' + file_read + '_qg.txt', delimiter=' ', skiprows=0, usecols = [0,2,3,4], unpack = True)
        hobs = 0.7
        ind = np.where(dn < 0)
        dn[ind] = 1e-9
        pDlog = np.log10(pD) +  3.0 * np.log10(hobs/h0)
        dnlog = np.log10(pD) - np.log10(dn)
        dulog = np.log10(du) - np.log10(pD)
        lm = lm -  2.0 * np.log10(hobs/h0)
        zappend.append((observation("Weaver+2023", lm, pDlog, abs(dulog), abs(dnlog), err_absolute=False), 'D'))

    # z0.5 obs
    z05obs = []
    add_weaver22_data(z05obs, file_read='0.2z0.5')
    add_weaver22_data(z05obs, file_read='0.5z0.8')

    # z1 obs
    z1obs = []
    add_weaver22_data(z1obs, file_read='0.8z1.1')
    add_weaver22_data(z1obs, file_read='1.1z1.5')


    #z2 obs
    z2obs = []
    add_weaver22_data(z2obs, file_read='2.0z2.5')
    add_weaver22_data(z2obs, file_read='1.5z2.0')

    # z3 obs
    z3obs = []
    add_weaver22_data(z3obs, file_read='3.0z3.5')
    add_weaver22_data(z3obs, file_read='2.5z3.0')


    # z4 obs
    z4obs = []
    add_weaver22_data(z4obs, file_read='3.5z4.5')

    # z5 obs
    z5obs = []
    add_weaver22_data(z5obs,file_read='4.5z5.5')


    return (z05obs, z1obs, z2obs, z3obs, z4obs, z5obs)

def load_smf_observations(obsdir, h0):
    # Thorne et al. (2021)
    def add_thorne21_data(zappend, file_read='z0.4'):
        lm, pD, dn, du = np.loadtxt(obsdir+'/mf/SMF/Thorne21/SMFvals_'+file_read+'.csv', delimiter=',', skiprows=1, usecols = [0,1,2,3], unpack = True)
        hobd = 0.7
        pDlog = np.log10(pD[::3]) +  3.0 * np.log10(hobs/h0)
        dnlog = np.log10(pD[::3]) - np.log10(dn[::3])
        dulog = np.log10(du[::3]) - np.log10(pD[::3])
        lm = lm[::3] -  2.0 * np.log10(hobs/h0)
        zappend.append((observation("Thorne+2021", lm, pDlog, abs(dulog), abs(dnlog), err_absolute=False), 's'))

    # Weaver et al. (2022; COSMOS2020)
    def add_weaver22_data(zappend, file_read='0.2z0.5', label = True):
        lm, pD, dn, du = np.loadtxt(obsdir+'/mf/SMF/COSMOS2020/SMF_Farmer_v2.1_' + file_read + '_total.txt', delimiter=' ', skiprows=0, usecols = [0,2,3,4], unpack = True)
        hobd = 0.7
        ind = np.where(dn < 0)
        dn[ind] = 1e-9
        pDlog = np.log10(pD) +  3.0 * np.log10(hobs/h0)
        dnlog = np.log10(pD) - np.log10(dn)
        dulog = np.log10(du) - np.log10(pD)
        lm = lm -  2.0 * np.log10(hobs/h0)
        zappend.append((observation("Weaver+2023" if label else None, lm, pDlog, abs(dulog), abs(dnlog), err_absolute=False), 'D'))


    # Driver al. (2022, z=0). Chabrier IMF
    z0obs = []
    z0obsPSO = []
    lm, p, dp = common.load_observation(obsdir, 'mf/SMF/GAMAIV_Driver22.dat', [0,1,2])
    hobs = 0.7
    xobs = lm + 2.0 * np.log10(hobs/h0)
    yobs = p - 3.0 * np.log10(hobs/h0)
    z0obs.append((observation("Driver+2022", xobs, yobs, dp, dp, err_absolute=False), 'o'))

    lm, p, dpdn, dpup = common.load_observation(obsdir, 'mf/SMF/SMF_Bernardi2013_SerExp.data', [0,1,2,3])
    xobs = lm + 2.0 * np.log10(hobs/h0)
    yobs = np.log10(p) - 3.0 * np.log10(hobs/h0)
    ydn = np.log10(p) - np.log10(p-dpdn)
    yup = np.log10(p+dpup) - np.log10(p) 
    z0obs.append((observation("Bernardi+2013", xobs, yobs, ydn, yup, err_absolute=False), 's'))


    lm, p, dpdn, dpup = common.load_observation(obsdir, 'mf/SMF/SMF_Li2009.dat', [0,1,2,3])
    xobs = lm - 2.0 * np.log10(hobs) + 2.0 * np.log10(hobs/h0)
    yobs = p + 3.0 * np.log10(hobs) - 3.0 * np.log10(hobs/h0)
    #z0obs.append((observation("Li&White+2009", xobs, yobs, abs(dpdn), dpup, err_absolute=False), 'd'))
    z0obsPSO.append((observation("Li+2009 (PSO)", xobs, yobs, abs(dpdn), dpup, err_absolute=False), 'd'))


    # Moustakas (Chabrier IMF), ['Moustakas+2013, several redshifts']
    zdnM13, lmM13, pM13, dp_dn_M13, dp_up_M13 = common.load_observation(obsdir, 'mf/SMF/SMF_Moustakas2013.dat', [0,3,5,6,7])
    xobsM13 = lmM13 + 2.0 * np.log10(hobs/h0)

    yobsM13 = np.full(xobsM13.shape, -999.) - 3.0 * np.log10(hobs/h0)
    lerrM13 = np.full(xobsM13.shape, -999.)
    herrM13 = np.full(xobsM13.shape, 999.)
    indx = np.where( pM13 < 1)
    yobsM13[indx] = (pM13[indx])
    indx = np.where( dp_dn_M13 > 0)
    lerrM13[indx]  = dp_dn_M13[indx] 
    indx = np.where( dp_up_M13 > 0)
    herrM13[indx]  = dp_up_M13[indx]

    # Muzzin (Kroupa IMF), ['Moustakas+2013, several redshifts']
    zdnMu13,zupMu13,lmMu13,pMu13,dp_dn_Mu13,dp_up_Mu13 = common.load_observation(obsdir, 'mf/SMF/SMF_Muzzin2013.dat', [0,1,2,4,5,5])
    # -0.09 corresponds to the IMF correction
    xobsMu13 = lmMu13 - 0.09 + 2.0 * np.log10(hobs/h0) 
    yobsMu13 = np.full(xobsMu13.shape, -999.) - 3.0 * np.log10(hobs/h0)
    lerrMu13 = np.full(xobsMu13.shape, -999.)
    herrMu13 = np.full(xobsMu13.shape, 999.)
    indx = np.where( pMu13 < 1)
    yobsMu13[indx] = (pMu13[indx])
    indx = np.where( dp_dn_Mu13 > 0)
    lerrMu13[indx]  = dp_dn_Mu13[indx] 
    indx = np.where( dp_up_Mu13 > 0)
    herrMu13[indx]  = dp_up_Mu13[indx]

    # z0.5 obs
    z05obs = []
    in_redshift = np.where(zdnMu13 == 0.5)
    z05obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z05obs, file_read='z0.51')
    add_weaver22_data(z05obs, file_read='0.2z0.5')
    add_weaver22_data(z05obs, file_read='0.5z0.8', label = False)


    # z1 obs
    z1obs = []
    in_redshift = np.where(zdnMu13 == 1)
    z1obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z1obs, file_read='z1.1')
    add_weaver22_data(z1obs, file_read='0.8z1.1', label = False)

    #z2 obs
    z2obs = []
    in_redshift = np.where(zupMu13 == 2.5)
    z2obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z2obs, file_read='z2')
    add_weaver22_data(z2obs, file_read='2.0z2.5')
    add_weaver22_data(z2obs, file_read='1.5z2.0', label = False)

    # z3 obs
    z3obs = []
    in_redshift = np.where(zupMu13 == 3.0)
    z3obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z3obs, file_read='z3')
    add_weaver22_data(z3obs, file_read='3.0z3.5')
    add_weaver22_data(z3obs, file_read='2.5z3.0', label = False)

    # z4 obs
    z4obs = []
    in_redshift = np.where(zupMu13 == 4.0)
    z4obs.append((observation("Muzzin+2013", xobsMu13[in_redshift], yobsMu13[in_redshift], lerrMu13[in_redshift], herrMu13[in_redshift], err_absolute=False), 'o'))
    add_thorne21_data(z4obs, file_read='z4')
    add_weaver22_data(z4obs, file_read='3.5z4.5')

    # z5 obs
    z5obs = []
    add_weaver22_data(z5obs,file_read='4.5z5.5')

    # z6 obs
    z6obs = []
    add_weaver22_data(z6obs,file_read='5.5z6.5')

    # z7 obs
    z7obs = []
    add_weaver22_data(z7obs,file_read='6.5z7.5')

    return (z0obs, z05obs, z1obs, z2obs, z3obs, z4obs, z5obs, z6obs, z7obs, z0obsPSO)

def plot_SMHM_z(plt, outdir, zlist, halo_mass_rel):

    def plot_moster13(ax, z, labels, label):
        #Moster et al. (2013) abundance matching SMHM relation
        M10 = 11.590
        M11 = 1.195
        N10 = 0.0351
        N11 = -0.0247
        beta10 = 1.376
        beta11 = -0.826
        gamma10 = 0.608
        gamma11 = 0.329
        M1 = pow(10.0, M10 + M11 * z/(z+1))
        N = N10 + N11 * z/(z+1)
        beta = beta10 + beta11 * z/(z+1)
        gamma = gamma10 + gamma11 * z/(z+1)

        mh = pow(10.0,xmf)
        m = mh * 2*N * pow (( pow(mh/M1, -beta ) + pow(mh/M1, gamma)), -1)

        if not labels:
            ax.plot(xmf,np.log10(m),'magenta', linestyle='dashed', linewidth=3)
        else:
            ax.plot(xmf,np.log10(m),'magenta', linestyle='dashed', linewidth=3, label=label)

    def plot_berhoozi13(ax, z, labels, label):
        a = 1.0/(1.0+z)
        nu = np.exp(-4*a*a)
        log_epsilon = -1.777 + (-0.006*(a-1)) * nu
        M1= 11.514 + ( - 1.793 * (a-1) - 0.251 * z) * nu
        alpha = -1.412 + 0.731 * nu * (a-1)
        delta = 3.508 + (2.608*(a-1)-0.043 * z) * nu
        gamma = 0.316 + (1.319*(a-1)+0.279 * z) * nu
        Min = xmf-M1
        fx = -np.log10(pow(10,alpha*Min)+1.0)+ delta * pow(np.log10(1+np.exp(Min)),gamma) / (1+np.exp(pow(10,-Min)))
        f = -0.3+ delta * pow(np.log10(2.0),gamma) / (1+np.exp(1))

        m = log_epsilon + M1 + fx - f

        if not labels:
            ax.plot(xmf,m, 'b', linestyle='dashdot',linewidth=3)
        else:
            ax.plot(xmf,m, 'b', linestyle='dashdot',linewidth=3, label=label)



    fig = plt.figure(figsize=(5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    xmin, xmax, ymin, ymax = 10.5, 14, 7, 13
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)


    subplots = (311, 312, 313)
    #zlist = (0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
    indices = [2,3,4]

    for i,idx in enumerate(indices):
        
        # z=0 ##################################
        ax = fig.add_subplot(subplots[i])
        if(i == 2):
            xtitle = xtit
        else:
            xtitle = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitle, ytit, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(zlist[idx]))

        ax.tick_params(labelsize=13)

        #Predicted SMHM
        ind = np.where(halo_mass_rel[idx,0,0,:] != 0)
        xplot = xmf[ind]
        yplot = halo_mass_rel[idx,0,0,ind]
        errdn = halo_mass_rel[idx,0,1,ind]
        errup = halo_mass_rel[idx,0,2,ind]
   
        ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='red', alpha = 0.5, interpolate=True)
        ax.errorbar(xplot, yplot[0], color='r', label='Passive galaxies (v2.0)')
        print("#SMHM relation of passive galaxies ")
        for a,b,c,d in zip(xplot, yplot[0],yplot[0]+errup[0],yplot[0]-errdn[0]):
            print(a,b,c,d, zlist[idx], 1)

        ind = np.where(halo_mass_rel[idx,2,0,:] != 0)
        xplot = xmf[ind]
        yplot = halo_mass_rel[idx,2,0,ind]
        errdn = halo_mass_rel[idx,2,1,ind]
        errup = halo_mass_rel[idx,2,2,ind]
   
        ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='blue', alpha = 0.5, interpolate=True)
        ax.errorbar(xplot, yplot[0], color='b', label = 'SF galaxies (v2.0)')
        print("#SMHM relation of SF galaxies ")
        for a,b,c,d in zip(xplot, yplot[0],yplot[0]+errup[0],yplot[0]-errdn[0]):
            print(a,b,c,d, zlist[idx], 0)
        ind = np.where(halo_mass_rel[idx,3,0,:] != 0)
        xplot = xmf[ind]
        yplot = halo_mass_rel[idx,3,0,ind]
        errdn = halo_mass_rel[idx,3,1,ind]
        errup = halo_mass_rel[idx,3,2,ind]
   
        ax.errorbar(xplot, yplot[0], color='k', label = 'all (v2.0)')

        #plot_moster13(ax, zlist[idx], True, 'Moster+13')
        #plot_berhoozi13(ax, zlist[idx], True, 'Behroozi+13')

        if(i == 0):
           common.prepare_legend(ax, ['r','b','k','magenta','b'], loc=2)

    plt.tight_layout()
    common.savefig(outdir, fig, 'SMHM_z_passivegals.pdf')


def plot_stellarmf_z(plt, outdir, obsdir, h0, hist_smf, hist_smf_err, hist_smf_30kpc, hist_smf_pass, hist_smf_pass_err):

    (z0obs, z05obs, z1obs, z2obs, z3obs, z4obs, z5obs, z6obs, z7obs, z0obsPSO) = load_smf_observations(obsdir, h0)
    
    PlotLagos18 = True
    def plot_lagos18_smf(ax, z):
        sm, z0, z0p5, z1, z2, z3, z4, z5, z6, z7 = common.load_observation(obsdir, 'Models/SharkVariations/SMF_Lagos18.dat', [0,1,2,3,4,5,6,7,8,9])
        if z == 0:
           y = z0
        elif z == 1:
           y = z0p5
        elif z == 2:
           y = z1
        elif z == 3:
           y = z2
        elif z == 4:
           y = z3
        elif z == 5:
           y = z4
        elif z == 6:
           y = z5
        elif z == 7:
           y = z6
        elif z == 8:
           y = z7
        ax.plot(sm, y, linestyle='dashed', color='black',label='Shark v1.1 (L18)' if z == 0 else None)

    
    fig = plt.figure(figsize=(11.7,11.7))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm M_{\\star}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 8, 13, -6, -1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (331, 332, 333, 334, 335, 336, 337, 338, 339)
    indeces = (0, 1, 2, 3, 4, 5, 6, 7, 8)
    zs = (0, 0.5, 1, 2, 3, 4, 5, 6, 7)
    observations = (z0obs, z05obs, z1obs, z2obs, z3obs, z4obs, z5obs, z6obs, z7obs)

    for subplot, idx, z, obs_and_markers in zip(subplots, indeces, zs, observations):

        ax = fig.add_subplot(subplot)
        if(idx ==0 or idx == 3 or idx == 6):
            ytitle = ytit
        else:
            ytitle = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(z))

        if( idx == 0):
            for obs, marker in z0obsPSO:
                common.errorbars(ax, obs.x, obs.y, obs.yerrdn, obs.yerrup, 'darkgreen',
                                   marker, err_absolute=obs.err_absolute, label=obs.label, markersize=6) 
        # Observations
        for obs, marker in obs_and_markers:
            common.errorbars(ax, obs.x, obs.y, obs.yerrdn, obs.yerrup, 'grey',
                                marker, err_absolute=obs.err_absolute, label=obs.label, markersize=4)

        # Predicted SMF
        y = hist_smf[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'r', label='Shark v2.0' if idx == 0 else None)
        #if idx == 0:
        #    y = hist_smf_err[idx,:]
        #    ind = np.where(y < 0.)
        #    ax.plot(xmf[ind],y[ind],'r', linestyle='dashdot', linewidth=2, label ='with error')

        if idx < 1:
            y = hist_smf_30kpc[idx,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'r', linestyle='dotted', linewidth=4, label ='Shark v2.0 (30kpc)'  if idx == 0 else None)
        if idx >= 1:
            y = hist_smf_err[idx,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'r', linestyle='dashdot', linewidth=3, label ='with 0.3dex error')
        plot_lagos18_smf(ax, idx)
       

        colors = []
        if idx == 0:
            colors = ['r','r','k']
        #elif idx == 1:
        #    colors += ['r']
        elif idx > 1:
            colors = ['r']
        if idx == 0:
           colors += ['darkgreen', 'grey', 'grey','grey']
        else:
           colors += ['grey', 'grey','grey']

        if idx ==0 or idx == 3:
           common.prepare_legend(ax, colors)

    plt.tight_layout()
    common.savefig(outdir, fig, 'stellarmf_z.pdf')

def plot_stellarmf_passive_z(plt, outdir, obsdir, h0, hist_smf, hist_smf_err, hist_smf_pass_cen, hist_smf_pass_sat):

    (z05obs, z1obs, z2obs, z3obs, z4obs, z5obs) = load_smf_passive_observations(obsdir, h0)
    PlotLagos18 = True
    def plot_lagos18_smf(ax, z, label=True):
        sm, z0, z0p5, z1, z2, z3, z4, z5 = common.load_observation(obsdir, 'Models/SharkVariations/SMF_Passive_Lagos18.dat', [0,1,2,3,4,5,6,7])
        y = z0
        if z == 1:
           y = z0p5
        elif z == 2:
           y = z1
        elif z == 3:
           y = z2
        elif z == 4:
           y = z3
        elif z == 5:
           y = z4
        elif z == 6:
           y = z5

        smc, z0c, z0p5c, z1c, z2c, z3c, z4c, z5c = common.load_observation(obsdir, 'Models/SharkVariations/SMF_Passive_Centrals_Lagos18.dat', [0,1,2,3,4,5,6,7])
        yc = z0c
        if z == 1:
           yc = z0p5c
        elif z == 2:
           yc = z1c
        elif z == 3:
           yc = z2c
        elif z == 4:
           yc = z3c
        elif z == 5:
           yc = z4c
        elif z == 6:
           yc = z5c
 
        sms, z0s, z0p5s, z1s, z2s, z3s, z4s, z5s = common.load_observation(obsdir, 'Models/SharkVariations/SMF_Passive_Satellites_Lagos18.dat', [0,1,2,3,4,5,6,7])
        ys = z0s
        if z == 1:
           ys = z0p5s
        elif z == 2:
           ys = z1s
        elif z == 3:
           ys = z2s
        elif z == 4:
           ys = z3s
        elif z == 5:
           ys = z4s
        elif z == 6:
           ys = z5s
        yall = np.log10(10**yc + 10**ys)
        ax.plot(smc, yall, linestyle='solid', color='black',label='Shark v1.1 (L18)' if label == True else None)
        ax.plot(smc, yc, linestyle='dashed', color='black',label='L18 cens with err' if label == True else None)
        ax.plot(sms, ys, linestyle='dotted', color='black',label='L18 sats with err' if label == True else None)


    fig = plt.figure(figsize=(11.7,8))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10}(\Phi/dlog_{10}{\\rm M_{\\star}}/{\\rm Mpc}^{-3} )$"
    xmin, xmax, ymin, ymax = 8, 13, -6.5, -1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (231, 232, 233, 234, 235, 236)
    indeces = (1, 2, 3, 4, 5, 6)
    zs = (0.5, 1, 2, 3, 4, 5)
    observations = (z05obs, z1obs, z2obs, z3obs, z4obs, z5obs)
    limits_obs = [7, 8.4, 9.4, 10, 10.1, 10.2, 10.45]

    for subplot, idx, z, obs_and_markers in zip(subplots, indeces, zs, observations):

        ax = fig.add_subplot(subplot)
        if(idx ==1 or idx == 4):
            ytitle = ytit
        else:
            ytitle = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle, locators=(0.1, 1, 0.1))
        ax.text(xleg, yleg, 'z=%s' % str(z))
        xlim = limits_obs[idx]
        ax.plot([xlim, xlim], [-6.5,-1], linestyle='dotted', color='grey')

        # Observations
        for obs, marker in obs_and_markers:
            if idx == 6:
                common.errorbars(ax, obs.x, obs.y, obs.yerrdn, obs.yerrup, 'grey',
                                marker, err_absolute=obs.err_absolute, label=obs.label, markersize=6)
            else:
                common.errorbars(ax, obs.x, obs.y, obs.yerrdn, obs.yerrup, 'grey',
                                marker, err_absolute=obs.err_absolute, markersize=6)


        # Predicted SMF
        y = hist_smf[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'r', label='Shark v2.0' if idx == 1 else None)
        y = hist_smf_err[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'r', linestyle='dashdot', linewidth=2, label ='0.3dex error' if idx == 1 else None)
        y = hist_smf_pass_cen[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'DarkOrange', label='centrals with err' if idx == 1 else None)
        y = hist_smf_pass_sat[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'YellowGreen', label='satellites with err' if idx == 1 else None)

        if(PlotLagos18):
            if(idx == 6):
               plot_lagos18_smf(ax, idx, label=True)
            else:
               plot_lagos18_smf(ax, idx, label=False)

        colors = []
        if idx == 1:
            colors = ['r','r','DarkOrange', 'YellowGreen', 'k']
        if idx == 6:
            colors=['k','k','k']
        colors += ['grey', 'grey','grey']

        if idx ==1:
           common.prepare_legend(ax, colors, loc = 3)
        elif idx == 6:
           common.prepare_legend(ax, colors, loc = 2)


    plt.tight_layout()
    common.savefig(outdir, fig, 'stellarmf_passive_z.pdf')


def prepare_data(hdf5_data, index, hist_smf, hist_smf_err, hist_smf_30kpc, hist_smf_pass, hist_smf_pass_err, hist_smf_pass_cen, hist_smf_pass_sat, sigma_ms, zlist, halo_mass_rel):

    (h0, volh, mdisk, mbulge, sfrd, sfrb, typeg, rstar_disk, rstar_bulge, mvir, mzd, mzb, mgd, mgb) = hdf5_data

    Zgas = ((mzd + mzb) / (mgd + mgb) / Zsun)

    bin_it = functools.partial(us.wmedians, xbins=xmf, nmin=10)

    ssfr = (sfrd + sfrb) / 1e9 / (mdisk + mbulge)
    ind = np.where(ssfr <= 1e-13)
    ssfr[ind] = 1e-13
    ssfr = np.log10(ssfr)

    mass          = np.zeros(shape = len(mdisk))
    ind = np.where((mdisk+mbulge) > 0.0)
    mass[ind] = np.log10(mdisk[ind] + mbulge[ind]) - np.log10(float(h0))
    H, _ = np.histogram(mass,bins=np.append(mbins,mupp))
    hist_smf[index,:] = hist_smf[index,:] + H
    ran_err = np.random.normal(0.0, 0.3, len(mass))
    mass_err = mass + ran_err
    H, _ = np.histogram(mass_err,bins=np.append(mbins,mupp))
    hist_smf_err[index,:] = hist_smf_err[index,:] + H

    ind = np.where(mass < 11)
    print("cosmic SFR at z", zlist[index], " is", np.log10(sum(sfrd[ind] + sfrb[ind]) / 1e9/h0/(volh / h0**3)))
 
    ind = np.where((mass > 0) & (ssfr <= -10.75))
    H, _ = np.histogram(mass[ind], bins=np.append(mbins,mupp))
    hist_smf_pass[index,:] = hist_smf_pass[index,:] + H

    ind = np.where((typeg > 0) & (mass > 10) & (ssfr <= -10.75) & (ssfr > -20) & (Zgas > 1e-10) & (Zgas < 1e10))
    print("Stats for satellites", np.log10(np.median(Zgas[ind])), np.median(sfrd[ind]/(sfrd[ind] + sfrb[ind])), " at redshift", zlist[index], len(sfrb[ind]))
    ind = np.where((typeg == 0) & (mass > 10) & (ssfr <= -10.75) & (ssfr > -20) & (Zgas > 1e-10) & (Zgas < 1e10))
    print("Stats for centrals", np.log10(np.median(Zgas[ind])), np.median(sfrd[ind]/(sfrd[ind] + sfrb[ind])), " at redshift", zlist[index], len(sfrb[ind]))

    #ind = np.where((mass > 0) & (ssfr <= -10.75) & (typeg == 0))
    #H, _ = np.histogram(mass[ind], bins=np.append(mbins,mupp))
    #hist_smf_pass_cen[index,:] = hist_smf_pass_cen[index,:] + H
    #ind = np.where((mass > 0) & (ssfr <= -10.75) & (typeg > 0))
    #H, _ = np.histogram(mass[ind], bins=np.append(mbins,mupp))
    #hist_smf_pass_sat[index,:] = hist_smf_pass_sat[index,:] + H

    if index == 0:
        scatter = 0.2
    else:
        scatter = 0.3
    ran_err = np.random.normal(0.0, scatter, len(mass))
    ssfr_err = ssfr + ran_err
    ind = np.where((mass_err > 0) & (ssfr_err <= -10.75))
    H, _ = np.histogram(mass_err[ind], bins=np.append(mbins,mupp))
    hist_smf_pass_err[index,:] = hist_smf_pass_err[index,:] + H

    ind = np.where((mass_err > 0) & (ssfr_err <= -10.75) & (typeg == 0))
    H, _ = np.histogram(mass_err[ind], bins=np.append(mbins,mupp))
    hist_smf_pass_cen[index,:] = hist_smf_pass_cen[index,:] + H
    halo_mass_rel[index,1,:] = bin_it(x=np.log10(mvir[ind]), y = mass_err[ind])

    ind = np.where((mass > 0) & (ssfr <= -10.75) & (typeg == 0))
    halo_mass_rel[index,0,:] = bin_it(x=np.log10(mvir[ind]), y = mass[ind])

    ind = np.where((mass_err > 0) & (ssfr_err <= -10.75) & (typeg > 0))
    H, _ = np.histogram(mass_err[ind], bins=np.append(mbins,mupp))
    hist_smf_pass_sat[index,:] = hist_smf_pass_sat[index,:] + H

    massd_30kpc = mdisk
    massb_30kpc = mbulge
    mass_30kpc = mdisk + mbulge
    #Calculate the stellar mass contained in 30pkpc, assuming an exponential profile for the disk and a Plummer profile for the bulge.
    ind = np.where((mdisk > 0.0)  & (rstar_disk > 0))
    massd_30kpc[ind] = mdisk[ind] * (1.0 - (1.0 + 30.0/(rstar_disk[ind]/1.67/h0 * MpcToKpc)) * np.exp(-30.0/(rstar_disk[ind]/1.67/h0 * MpcToKpc)))
    ind = np.where((mbulge > 0.0)  & (rstar_bulge > 0))
    massb_30kpc[ind] = mbulge[ind] * pow(30.0, 3.0) / pow((pow(30.0, 2.0) + pow(rstar_bulge[ind]/1.3/h0 * MpcToKpc, 2.0)), 3.0/2.0)

    ind = np.where((massd_30kpc + massb_30kpc) > 0)
    mass_30kpc[ind] = np.log10(massd_30kpc[ind] + massb_30kpc[ind]) - np.log10(float(h0))
    H, _ = np.histogram(mass_30kpc,bins=np.append(mbins,mupp))
    hist_smf_30kpc[index,:] = hist_smf_30kpc[index,:] + H

    if volh > 0:
        vol = volh/pow(h0,3.)  # In Mpc^3
        hist_smf[index,:]  = hist_smf[index,:]/vol/dm
        hist_smf_err[index,:]  = hist_smf_err[index,:]/vol/dm
        hist_smf_30kpc[index,:]  = hist_smf_30kpc[index,:]/vol/dm
        hist_smf_pass[index,:]  = hist_smf_pass[index,:]/vol/dm
        hist_smf_pass_err[index,:]  = hist_smf_pass_err[index,:]/vol/dm
        hist_smf_pass_cen[index,:]  = hist_smf_pass_cen[index,:]/vol/dm
        hist_smf_pass_sat[index,:]  = hist_smf_pass_sat[index,:]/vol/dm


    sfrt = (sfrd+ sfrb)/1e9/h0
    mt = (mdisk + mbulge)/h0


    bin_it = functools.partial(us.wmedians, xbins=xmf)
    bin_it_2 = functools.partial(us.wmedians, xbins=xmf2)

    #select main sequence galaxies
    ind = np.where((sfrt > 0) & (mt >= 7e8) & (mt <= 1e10) & (typeg == 0))
    sfrin = np.log10(sfrt[ind])
    mtin = np.log10(mt[ind])

    ms_med = bin_it(x=mtin, y=sfrin)
    pos = np.where(ms_med[0,:] != 0)
    yin = ms_med[0,pos]
    ms_fit = np.polyfit(xmf[pos], yin[0], 2)
   
    dist_ms = np.log10(sfrt) - (ms_fit[0] * np.log10(mt)**2 + ms_fit[1] * np.log10(mt) + ms_fit[2])

    for i,m in enumerate(xmf2):
        ind = np.where((sfrt > 0) & (mt >= 10**(m-dm2/2.0)) & (mt < 10**(m+dm2/2.0))  & (dist_ms > -0.75) & (dist_ms < 0.75))
        sigma_ms[index,i] = np.std(np.log10(sfrt[ind]))
   
    
    ind = np.where((mass > 0) & (typeg == 0))
    halo_mass_rel[index,3,:] = bin_it(x=np.log10(mvir[ind]), y = mass[ind])

    ind = np.where((mass > 0) & (dist_ms > -0.3) & (typeg == 0))
    halo_mass_rel[index,2,:] = bin_it(x=np.log10(mvir[ind]), y = mass[ind])

    #print("#Main sequence scatter at z:", zlist[index])
    #print("#log10(Mstar/Msun) std(log10(SFR))")

    #for a,b in zip(xmf2, sigma_ms[index,:]):
    #    print(a,b)
 

    return mass

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist = (0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0) #(0.15, 0.25, 0.4, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875)

    plt = common.load_matplotlib()

    # Histograms
    hist_smf       = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_err   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_30kpc  = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_pass   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_pass_cen   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_pass_sat   = np.zeros(shape = (len(zlist), len(mbins)))
    hist_smf_pass_err = np.zeros(shape = (len(zlist), len(mbins)))
    sigma_ms       = np.zeros(shape = (len(zlist), len(xmf2)))
    halo_mass_rel  = np.zeros(shape = (len(zlist), 4, 3, len(mbins)))

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'sfr_disk', 'sfr_burst', 'type', 'rstar_disk', 'rstar_bulge', 'mvir_hosthalo',
                           'mgas_metals_disk', 'mgas_metals_bulge', 'mgas_disk', 'mgas_bulge')}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        mass = prepare_data(hdf5_data, index, hist_smf, hist_smf_err, hist_smf_30kpc, hist_smf_pass, hist_smf_pass_err, hist_smf_pass_cen, hist_smf_pass_sat, sigma_ms, zlist, halo_mass_rel)
        h0 = hdf5_data[0]

    # Take logs
    def take_log(array):
        ind = np.where(array > 0.)
        array[ind] = np.log10(array[ind])
    take_log(hist_smf)
    take_log(hist_smf_err)
    take_log(hist_smf_30kpc)
    take_log(hist_smf_pass)
    take_log(hist_smf_pass_err)
    take_log(hist_smf_pass_cen)
    take_log(hist_smf_pass_sat)

    plot_stellarmf_z(plt, outdir, obsdir, h0, hist_smf, hist_smf_err, hist_smf_30kpc, hist_smf_pass, hist_smf_pass_err)
    plot_stellarmf_passive_z(plt, outdir, obsdir, h0, hist_smf_pass, hist_smf_pass_err, hist_smf_pass_cen, hist_smf_pass_sat)
    plot_SMHM_z(plt, outdir, zlist, halo_mass_rel)
    #print("#SMF passive galaxies")
    #for a,b,c,d,e,f,g,h in zip(xmf, hist_smf_pass_sat[0,:], hist_smf_pass_sat[1,:], hist_smf_pass_sat[2,:], hist_smf_pass_sat[3,:], hist_smf_pass_sat[4,:], hist_smf_pass_sat[5,:], hist_smf_pass_sat[6,:]):
    #    print(a,b,c,d,e,f,g,h)

if __name__ == '__main__':
    main(*common.parse_args())
