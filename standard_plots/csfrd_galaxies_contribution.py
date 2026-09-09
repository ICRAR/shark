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
XH = 0.72
MpcToKpc = 1e3
Lsun = 3.828e-7 #in 1e40 erg/s
conv_ir_sfr = 10**3.41  #conversion from 1e40erg/s to Msun/yr (Kennicutt 2013)
conv_uv_sfr = 10**3.35  #conversion from 1e40erg/s to Msun/yr (Kennicutt 2013)
dist_cm2_pc2 = 4.0 * np.pi * (10 * 3.086e18)**2 #distance to 10 parsec [in cm] for luminosity conversion
uv_1500_to_Hz = 1.998e15 #conversion from 1500 Angstrom to Hz

ab_to_lsun = 10**(np.log10(dist_cm2_pc2) + np.log10(uv_1500_to_Hz) - 40.0)

xmf_edges = np.array([9, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 12.5])
xmf = np.array([9.125, 9.375, 9.625, 9.875, 10.125, 10.375, 10.625, 10.875, 11.125, 11.375, 11.675])

zsun = 0.0189
#model of Mattson et al. (2014) for the dependence of the dust-to-metal mass ratio and metallicity X/H.
corrfactor_dm = 2.0
polyfit_dm = [ 0.00544948, 0.00356938, -0.07893235,  0.05204814,  0.49353238]

#choose dust model between mm14, rr14 and constdust
m14 = False
rr14 = True
constdust = False
rr14xcoc = False

# compute dust masses
def dust_mass(mz, mg, h0):
    md = np.zeros(shape = len(mz))
    ind = np.where((mz > 0) & (mg > 0))
    XHd = np.log10(mz[ind]/mg[ind]/zsun)
    if(m14 == True):
        DToM = (polyfit_dm[0] * XHd**4.0 + polyfit_dm[1] * XHd**3.0 + polyfit_dm[2] * XHd**2.0 + polyfit_dm[3] * XHd + polyfit_dm[4])/corrfactor_dm
        DToM = np.clip(DToM, 1e-6, 0.5)
        md[ind] = mz[ind]/h0 * DToM
        DToM_MW = polyfit_dm[4]/corrfactor_dm
    elif(rr14 == True):
         y = np.zeros(shape = len(XHd))
         highm = np.where(XHd > -0.59)
         y[highm] = 10.0**(2.21 - XHd[highm]) #gas-to-dust mass ratio
         lowm = np.where(XHd <= -0.59)
         y[lowm] = 10.0**(0.96 - (3.1) * XHd[lowm]) #gas-to-dust mass ratio
         DToM = 1.0 / y / (mz[ind]/mg[ind])
         DToM = np.clip(DToM, 1e-6, 1)
         md[ind] = mz[ind]/h0 * DToM
         DToM_MW = 1.0 / (10.0**(2.21)) / zsun
    elif(rr14xcoc == True):
         y = np.zeros(shape = len(XHd))
         highm = np.where(XHd > -0.15999999999999998)
         y[highm] = 10.0**(2.21 - XHd[highm]) #gas-to-dust mass ratio
         lowm = np.where(XHd <= -0.15999999999999998)
         y[lowm] = 10.0**(1.66 - 4.43 * XHd[lowm]) #gas-to-dust mass ratio
         DToM = 1.0 / y / (mz[ind]/mg[ind])
         DToM = np.clip(DToM, 1e-6, 1)
         md[ind] = mz[ind]/h0 * DToM
         DToM_MW = 1.0 / (10.0**(2.21)) / zsun
    elif(constdust == True):
         md[ind] = 0.33 * mz[ind]/h0
         DToM_MW = 0.33

    return (md, DToM_MW)

def plot_obs_zalava21(ax, h0):

    ### IR data
    hobs = 0.7
    #constants are to convert IMF to a Chabrier IMF
    zM14, zlM14, zuM14, lgrhoM14, lgrhoM14_up, lgrhoM14dn = np.loadtxt('../data/Global/sfrd_literature/madau14a_irdata.txt', unpack = True, usecols=[0,1,2,3,4,5])
    ax.errorbar(zM14,lgrhoM14 + np.log10(0.63) + np.log10(hobs/h0),xerr=[zM14-zlM14, zuM14-zM14], yerr=[lgrhoM14dn,lgrhoM14_up], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5, label='Zavala+21 (IR)')
    zC14, lgrhoC14, lgrhoC14_up, lgrhoC14dn = np.loadtxt('../data/Global/sfrd_literature/casey14a_irdata.txt', unpack = True, usecols=[0,1,2,3])
    #ax.errorbar(zC14,lgrhoC14 + np.log10(0.55) + np.log10(hobs/h0),yerr=[abs(lgrhoC14-lgrhoC14dn),abs(lgrhoC14_up-lgrhoC14)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zC21, zlC21, zuC21, lgrhoC21, lgrhoC21_dn, lgrhoC21_up = np.loadtxt('../data/Global/sfrd_literature/casey2021.txt', unpack = True, usecols=[0,1,2,3,4,5])
    lgrhoC21 = np.log10(lgrhoC21)
    lgrhoC21_dn = np.log10(lgrhoC21_dn)
    lgrhoC21_up = np.log10(lgrhoC21_up)
    ax.errorbar(zC21,lgrhoC21 + np.log10(hobs/h0), xerr=[abs(zC21 - zlC21), abs(zC21 - zuC21)], yerr=[abs(lgrhoC21-lgrhoC21_dn),abs(lgrhoC21_up-lgrhoC21)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zC13, lgrhoC13, lgrhoC13_dn, lgrhoC13_up = np.loadtxt('../data/Global/sfrd_literature/casey2013.txt', unpack = True, usecols=[0,3,4,5])
    ax.errorbar(zC13,lgrhoC13 + np.log10(0.55) + np.log10(hobs/h0),yerr=[abs(lgrhoC13-lgrhoC13_dn),abs(lgrhoC13_up-lgrhoC13)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zD20, lgrhoD20, lgrhoD20_dn, lgrhoD20_up = np.loadtxt('../data/Global/sfrd_literature/dudzeviciute2020.txt', unpack = True, usecols=[0,3,4,5])
    lgrhoD20 = np.log10(lgrhoD20)
    lgrhoD20_dn = np.log10(lgrhoD20_dn)
    lgrhoD20_up = np.log10(lgrhoD20_up)
    ax.errorbar(zD20,lgrhoD20 + np.log10(hobs/h0), yerr=[abs(lgrhoD20-lgrhoD20_dn),abs(lgrhoD20_up-lgrhoD20)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zG20, zlG20, zuG20, lgrhoG20, lgrhoG20_dn, lgrhoG20_up = np.loadtxt('../data/Global/sfrd_literature/gruppioni2020.txt', unpack = True, usecols=[0,1,2,3,4,5])
    lgrhoG20 = np.log10(lgrhoG20)
    lgrhoG20_dn = np.log10(lgrhoG20_dn)
    lgrhoG20_up = np.log10(lgrhoG20_up)
    ax.errorbar(zG20,lgrhoG20 + np.log10(hobs/h0), xerr=[abs(zG20 - zlG20), abs(zG20 - zuG20)], yerr=[abs(lgrhoG20-lgrhoG20_dn),abs(lgrhoG20_up-lgrhoG20)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zL20, zlL20, zuL20, lgrhoL20, lgrhoL20_dn, lgrhoL20_up = np.loadtxt('../data/Global/sfrd_literature/lim2020.txt', unpack = True, usecols=[0,1,2,3,4,5])
    lgrhoL20 = np.log10(lgrhoL20)
    lgrhoL20_dn = np.log10(lgrhoL20_dn)
    lgrhoL20_up = np.log10(lgrhoL20_up)
    ax.errorbar(zL20,lgrhoL20 + np.log10(hobs/h0), xerr=[abs(zL20 - zlL20), abs(zL20 - zuL20)], yerr=[abs(lgrhoL20-lgrhoL20_dn),abs(lgrhoL20_up-lgrhoL20)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zL18, zlL18, zuL18, lgrhoL18, lgrhoL18_dn, lgrhoL18_up = np.loadtxt('../data/Global/sfrd_literature/liu2018_err.txt', unpack = True, usecols=[0,1,2,3,4,5])
    lgrhoL18 = np.log10(lgrhoL18)
    lgrhoL18_dn = np.log10(lgrhoL18_dn)
    lgrhoL18_up = np.log10(lgrhoL18_up)
    ax.errorbar(zL18,lgrhoL18 + np.log10(hobs/h0), xerr=[abs(zL18 - zlL18), abs(zL18 - zuL18)], yerr=[abs(lgrhoL18-lgrhoL18_dn),abs(lgrhoL18_up-lgrhoL18)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zS14, zlS14, zuS14, lgrhoS14, lgrhoS14_dn, lgrhoS14_up = np.loadtxt('../data/Global/sfrd_literature/swinbank2014.txt', unpack = True, usecols=[0,1,2,3,4,5])
    lgrhoS14 = np.log10(lgrhoS14)
    lgrhoS14_dn = np.log10(lgrhoS14_dn)
    lgrhoS14_up = np.log10(lgrhoS14_up)
    ax.errorbar(zS14,lgrhoS14 + np.log10(hobs/h0), xerr=[abs(zS14 - zlS14), abs(zS14 - zuS14)], yerr=[abs(lgrhoS14-lgrhoS14_dn),abs(lgrhoS14_up-lgrhoS14)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zW19, zlW19, zuW19, lgrhoW19, lgrhoW19_dn, lgrhoW19_up = np.loadtxt('../data/Global/sfrd_literature/wangL2019.txt', unpack = True, usecols=[0,1,2,3,4,5])
    lgrhoW19 = np.log10(lgrhoW19)
    lgrhoW19_dn = np.log10(lgrhoW19_dn)
    lgrhoW19_up = np.log10(lgrhoW19_up)
    ax.errorbar(zW19,lgrhoW19 + np.log10(hobs/h0), xerr=[abs(zW19 - zlW19), abs(zW19 - zuW19)], yerr=[abs(lgrhoW19-lgrhoW19_dn),abs(lgrhoW19_up-lgrhoW19)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zK17, zlK17, zuK17, lgrhoK17, lgrhoK17_dn, lgrhoK17_up = np.loadtxt('../data/Global/sfrd_literature/koprowski2017.txt', unpack = True, usecols=[0,1,2,3,4,5])
    lgrhoK17 = np.log10(lgrhoK17)
    lgrhoK17_dn = np.log10(lgrhoK17_dn)
    lgrhoK17_up = np.log10(lgrhoK17_up)
    ax.errorbar(zK17,lgrhoK17 + np.log10(hobs/h0), xerr=[abs(zK17 - zlK17), abs(zK17 - zuK17)], yerr=[abs(lgrhoK17-lgrhoK17_dn),abs(lgrhoK17_up-lgrhoK17)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zB17, zlB17, zuB17, lgrhoB17, lgrhoB17_dn, lgrhoB17_up = np.loadtxt('../data/Global/sfrd_literature/bourne2017.txt', unpack = True, usecols=[0,1,2,3,4,5])
    lgrhoB17 = np.log10(lgrhoB17)
    lgrhoB17_dn = np.log10(lgrhoB17_dn)
    lgrhoB17_up = np.log10(lgrhoB17_up)
    ax.errorbar(zB17,lgrhoB17 + np.log10(hobs/h0), xerr=[abs(zB17 - zlB17), abs(zB17 - zuB17)], yerr=[abs(lgrhoB17-lgrhoB17_dn),abs(lgrhoB17_up-lgrhoB17)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zZ21, zlZ21, zuZ21, lgrhoZ21, lgrhoZ21_dn, lgrhoZ21_up = np.loadtxt('../data/Global/sfrd_literature/zavala2021.txt', unpack = True, usecols=[0,1,2,3,4,5])
    lgrhoZ21 = np.log10(lgrhoZ21)
    lgrhoZ21_dn = np.log10(lgrhoZ21_dn)
    lgrhoZ21_up = np.log10(lgrhoZ21_up)
    ax.errorbar(zZ21,lgrhoZ21 + np.log10(hobs/h0), xerr=[abs(zZ21 - zlZ21), abs(zZ21 - zuZ21)], yerr=[abs(lgrhoZ21-lgrhoZ21_dn),abs(lgrhoZ21_up-lgrhoZ21)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zB20, zlB20, zuB20, lgrhoB20, lgrhoB20_dn, lgrhoB20_up = np.loadtxt('../data/Global/sfrd_literature/bouwens2020.txt', unpack = True, usecols=[0,1,2,3,4,5])
    lgrhoB20 = np.log10(lgrhoB20)
    lgrhoB20_dn = np.log10(lgrhoB20_dn)
    lgrhoB20_up = np.log10(lgrhoB20_up)
    ax.errorbar(zB20,lgrhoB20 + np.log10(hobs/h0), xerr=[abs(zB20 - zlB20), abs(zB20 - zuB20)], yerr=[abs(lgrhoB20-lgrhoB20_dn),abs(lgrhoB20_up-lgrhoB20)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)
    zA22, zlA22, zuA22, lgrhoA22, lgrhoA22_dn, lgrhoA22_up = np.loadtxt('../data/Global/sfrd_literature/algera2022.txt', unpack = True, usecols=[0,1,2,3,4,5])
    lgrhoA22 = np.log10(lgrhoA22)
    lgrhoA22_dn = np.log10(lgrhoA22_dn)
    lgrhoA22_up = np.log10(lgrhoA22_up)
    ax.errorbar(zA22,lgrhoA22 + np.log10(hobs/h0), xerr=[abs(zA22 - zlA22), abs(zA22 - zuA22)], yerr=[abs(lgrhoA22-lgrhoA22_dn),abs(lgrhoA22_up-lgrhoA22)], ls='None', mfc='darkred', ecolor = 'darkred', mec='darkred',marker='s', alpha = 0.5)

    zZ21, lgrhoZ21, lgrhoZ21_dn, lgrhoZ21_up = np.loadtxt('../data/Global/sfrd_literature/zavala2021_uv_compilation.txt', unpack = True, usecols=[0,1,2,3])
    ax.errorbar(zZ21,lgrhoZ21 + np.log10(hobs/h0), yerr=[abs(lgrhoZ21-lgrhoZ21_dn),abs(lgrhoZ21_up-lgrhoZ21)], ls='None', mfc='blue', ecolor = 'blue', mec='blue',marker='o', alpha = 0.5, label='Zavala+21 (UV)')
    zL18, zlL18, zuL18, lgrhoL18, lgrhoL18_dn = np.loadtxt('../data/Global/SFRD_UV_Chemerynska.dat', unpack = True, usecols=[0,1,2,3,4])
    ax.errorbar(zL18,lgrhoL18 + np.log10(hobs/h0), xerr=[abs(zL18 - zlL18), abs(zL18 - zuL18)], yerr=[lgrhoL18_dn, lgrhoL18_dn], ls='None', mfc='blue', ecolor = 'blue', mec='blue',marker='o', alpha = 0.5)


def plot_SFRD(plt, outdir, obsdir, csfrd, zlist, h0):


    #comparison with micro SURFS
    fig = plt.figure(figsize=(14,5.5))

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(CSFRD/ M_{\odot}\,yr^{-1}\,cMpc^{-3})$"

    ax = fig.add_subplot(131)
    plt.subplots_adjust(left=0.15)

    common.prepare_ax(ax, 0, 15, -6, -0.3, xtit, ytit, locators=(1, 1, 1, 1))
    #plt.xscale('log')
    ax.text(5,-0.7,"Total", fontsize=14)
    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    reddnM14, redupM14, sfrM14, sfrM14errup, sfrM14errdn = common.load_observation(obsdir, 'Global/SFRD_Madau14.dat', [0,1,2,3,4])
    #authors assume a Salpeter IMF, so a correction of np.log10(0.63) is necessary.
    sfrM14errdn = abs(sfrM14errdn)
    hobs = 0.7
    sfrM14 = sfrM14 + np.log10(pow(hobs/h0, 2.0)) + np.log10(0.63)
    #ax.errorbar((reddnM14 + redupM14) / 2.0, sfrM14, xerr=abs(redupM14-reddnM14)/2.0, yerr=[sfrM14errdn, sfrM14errup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D', markersize=1.5, label='Madau+14')

    #D'Silva+23 (Chabrier IMF)
    sfrD23, sfrD23errup, sfrD23errdn, zD23, zD23errup, zD23errdn = common.load_observation(obsdir, 'Global/DSilva23_sfr.dat', [0,1,2,3,4,5])
    hobs = 0.7
    yobsD23 = sfrD23 + np.log10(hobs/h0)

    #Adams23 (Chabrier IMF)
    zA23, sfrA23, sfrA23errdn, sfrA23errup = common.load_observation(obsdir, 'Global/Adams23_CSFRDCompilation.dat', [0,1,2,3])
    sfrA23errdn = sfrA23 - sfrA23errdn #make them relative errors
    sfrA23errup = sfrA23errup - sfrA23
    hobs = 0.7
    yobsA23 = sfrA23 + np.log10(hobs/h0)
    ax.errorbar(zA23, yobsA23, yerr=[sfrA23errdn, sfrA23errup], ls='None', mfc='None', ecolor = 'darkgreen', mec='darkgreen',marker='s', label='Adams+24')

    #D'Silva+25 (Chabrier IMF)
    zD25, zD25dn, zD25up, sfrD25, sfrD25err = common.load_observation(obsdir, 'Global/DSilva25_sfr.dat', [0,1,2,3,4])
    zD25errdn = zD25-zD25dn
    zD25errup = zD25up - zD25

    hobs = 0.7
    yobsD25 = sfrD25 + np.log10(hobs/h0)

    ax.errorbar(zD23, yobsD23, xerr=[zD23errdn, zD23errup], yerr=[sfrD23errdn, sfrD23errup], ls='None', mfc='None', ecolor = 'navy', mec='navy',marker='o', label ='D\'Silva+23,25')
    ax.errorbar(zD25, yobsD25, xerr=[zD25errdn, zD25errup], yerr=sfrD25err, ls='None', mfc='None', ecolor = 'navy', mec='navy',marker='o')

    #Bouwens 26 (Chabrier IMF)
    zB26, sfrB26, sfrB26errdn, sfrB26errup = common.load_observation(obsdir, 'Global/SFRD_Compilation_Bouwens26.dat', [0,1,2,3])
    sfrB26errdn = sfrB26 - sfrB26errdn #make them relative errors
    sfrB26errup = sfrB26errup - sfrB26
    hobs = 0.7
    yobsB26 = sfrB26 + np.log10(hobs/h0)
    #ax.errorbar(zB26, yobsB26, yerr=[sfrB26errdn, sfrB26errup], ls='None', mfc='None', ecolor = 'darkgreen', mec='darkgreen',marker='s',label='Bouwens26 (comp)')

    #Driver (Chabrier IMF)
    redD17d, redD17u, sfrD17, err1, err2, err3, err4 = common.load_observation(obsdir, 'Global/Driver18_sfr.dat', [0,1,2,3,4,5,6])
    hobs = 0.7
    xobsD17 = (redD17d+redD17u)/2.0
    yobsD17 = sfrD17 + np.log10(hobs/h0)
    errD17 = yobsD17*0. - 999.
    errD17 = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))
    ax.errorbar(xobsD17, yobsD17, yerr=[errD17,errD17], ls='None', mfc='None', ecolor = 'darkorange', mec='darkorange',marker='o', label = 'Driver+18')

    ax.plot(zlist, np.log10(csfrd[:,0,0]), 'k', linewidth=3, label ='intr.')
    ax.plot(zlist, np.log10(csfrd[:,0,1]), 'k', linewidth=1, label ='eith err.')

    common.prepare_legend(ax, ['k','purple','olive','k','k','k'], loc=3) #bbox_to_anchor=(0.52, 0.47))

    ax = fig.add_subplot(132)
    plt.subplots_adjust(left=0.15)
    ax.text(3,-0.7,"Star formation in disks", fontsize=14)

    common.prepare_ax(ax, 0, 15, -6, -0.3, xtit, " ", locators=(1, 1, 1, 1))
    #plt.xscale('log')

    ax.errorbar(zD23, yobsD23, xerr=[zD23errdn, zD23errup], yerr=[sfrD23errdn, sfrD23errup], ls='None', mfc='None', ecolor = 'navy', mec='navy',marker='o', label ='D\'Silva+23,25')
    ax.errorbar(zD25, yobsD25, xerr=[zD25errdn, zD25errup], yerr=sfrD25err, ls='None', mfc='None', ecolor = 'navy', mec='navy',marker='o')

    ax.errorbar(zA23, yobsA23, yerr=[sfrA23errdn, sfrA23errup], ls='None', mfc='None', ecolor = 'darkgreen', mec='darkgreen',marker='s', label = 'Adams+24')
    ax.errorbar(xobsD17, yobsD17, yerr=[errD17,errD17], ls='None', mfc='None', ecolor = 'darkorange', mec='darkorange',marker='o', label = 'Driver+18')
    #ax.errorbar(zB26, yobsB26, yerr=[sfrB26errdn, sfrB26errup], ls='None', mfc='None', ecolor = 'darkgreen', mec='darkgreen',marker='s',label='Bouwens26 (comp)')

    ax.plot(zlist, np.log10(csfrd[:,1,0]), 'k', linewidth=3, label ='intr.')
    ax.plot(zlist, np.log10(csfrd[:,1,1]), 'k', linewidth=1, label ='eith err.')

    common.prepare_legend(ax, ['k','purple','olive', 'k','k','k'], loc=3) #bbox_to_anchor=(0.52, 0.47))

    ax = fig.add_subplot(133)
    plt.subplots_adjust(left=0.15)
    ax.text(5,-0.7,"Starbursts", fontsize=14)

    common.prepare_ax(ax, 0, 15, -6, -0.3, xtit, " ", locators=(1, 1, 1, 1))
    #plt.xscale('log')

    ax.errorbar(zD23, yobsD23, xerr=[zD23errdn, zD23errup], yerr=[sfrD23errdn, sfrD23errup], ls='None', mfc='None', ecolor = 'navy', mec='navy',marker='o', label ='D\'Silva+23,25')
    ax.errorbar(zD25, yobsD25, xerr=[zD25errdn, zD25errup], yerr=sfrD25err, ls='None', mfc='None', ecolor = 'navy', mec='navy',marker='o')

    ax.errorbar(zA23, yobsA23, yerr=[sfrA23errdn, sfrA23errup], ls='None', mfc='None', ecolor = 'darkgreen', mec='darkgreen',marker='s', label = 'Adams+24')
    ax.errorbar(xobsD17, yobsD17, yerr=[errD17,errD17], ls='None', mfc='None', ecolor = 'darkorange', mec='darkorange',marker='o', label = 'Driver+18')
    #ax.errorbar(zB26, yobsB26, yerr=[sfrB26errdn, sfrB26errup], ls='None', mfc='None', ecolor = 'darkgreen', mec='darkgreen',marker='s',label='Bouwens26 (comp)')

    ax.plot(zlist, np.log10(csfrd[:,2,0]), 'k', linewidth=3, label ='intr.')
    ax.plot(zlist, np.log10(csfrd[:,2,1]), 'k', linewidth=1, label ='eith err.')

    common.prepare_legend(ax, ['k','purple','olive', 'k','k','k'], loc=3) #bbox_to_anchor=(0.52, 0.47))
    plt.tight_layout()
    common.savefig(outdir, fig, "cosmic_sfr_compmicro_galcalculated.pdf")

def plot_SFRD_from_lum(plt, outdir, obsdir, csfrd, zlist, h0, file_name):

    #comparison with micro SURFS
    fig = plt.figure(figsize=(14,5.5))

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(CSFRD/ M_{\odot}\,yr^{-1}\,cMpc^{-3})$"

    ax = fig.add_subplot(131)
    plt.subplots_adjust(left=0.15)

    common.prepare_ax(ax, 0, 15, -6, -0.3, xtit, ytit, locators=(1, 1, 1, 1))
    #plt.xscale('log')
    ax.text(5,-0.7,"Total (L800)", fontsize=14)

    z_rr14steep, gcsfr_un, gcsfr_at, dcsfr_un, dcsfr_at, scsfr_un, scsfr_at = np.loadtxt('../data/Models/SharkVariations/CSFRD_Shark2_L800_eagle-rr14-steep.txt', unpack=True, usecols=[0,1,2,3,4,5,6])

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    #reddnM14, redupM14, sfrM14, sfrM14errup, sfrM14errdn = common.load_observation(obsdir, 'Global/SFRD_Madau14.dat', [0,1,2,3,4])
    ##authors assume a Salpeter IMF, so a correction of np.log10(0.63) is necessary.
    #sfrM14errdn = abs(sfrM14errdn)
    #hobs = 0.7
    #sfrM14 = sfrM14 + np.log10(pow(hobs/h0, 2.0)) + np.log10(0.63)
    #ax.errorbar((reddnM14 + redupM14) / 2.0, sfrM14, xerr=abs(redupM14-reddnM14)/2.0, yerr=[sfrM14errdn, sfrM14errup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D', markersize=1.5, label='Madau+14')

    #D'Silva+23 (Chabrier IMF)
    #sfrD23, sfrD23errup, sfrD23errdn, zD23, zD23errup, zD23errdn = common.load_observation(obsdir, 'Global/DSilva23_sfr.dat', [0,1,2,3,4,5])
    #hobs = 0.7
    #yobsD23 = sfrD23 + np.log10(hobs/h0)

    #Adams23 (Chabrier IMF)
    zA23, sfrA23, sfrA23errdn, sfrA23errup = common.load_observation(obsdir, 'Global/Adams23_CSFRDCompilation.dat', [0,1,2,3])
    sfrA23errdn = sfrA23 - sfrA23errdn #make them relative errors
    sfrA23errup = sfrA23errup - sfrA23
    hobs = 0.7
    yobsA23 = sfrA23 + np.log10(hobs/h0)
    ax.errorbar(zA23, yobsA23, yerr=[sfrA23errdn, sfrA23errup], ls='None', mfc='darkgreen', ecolor = 'darkgreen', mec='darkgreen',marker='*', label='Adams+24 (JWST)')

    #Zavala+21
    uv_zZ21, uv_sfrZ21 = common.load_observation(obsdir, 'Global/SFRD_Compilation_UV_Zavala21.dat', [0,1])
    hobs = 0.7
    uv_yobsZ21 = uv_sfrZ21 + np.log10(hobs/h0)
    ir_zZ21, ir_sfrZ21 = common.load_observation(obsdir, 'Global/SFRD_Compilation_IR_Zavala21.dat', [0,1])
    hobs = 0.7
    ir_yobsZ21 = ir_sfrZ21 + np.log10(hobs/h0)
    #ax.plot(uv_zZ21, uv_yobsZ21, color='blue', marker='s', ls='None', alpha = 0.5, label = 'Zavala+21 (UV)')
    #ax.plot(ir_zZ21, ir_yobsZ21, color='darkred', marker='o', ls='None', alpha= 0.5, label = 'Zavala+21 (IR)')

    plot_obs_zalava21(ax, h0)
    #ax.errorbar(zD23, yobsD23, xerr=[zD23errdn, zD23errup], yerr=[sfrD23errdn, sfrD23errup], ls='None', mfc='None', ecolor = 'navy', mec='navy',marker='o', label ='D\'Silva+23,25')
    #ax.errorbar(zD25, yobsD25, xerr=[zD25errdn, zD25errup], yerr=sfrD25err, ls='None', mfc='None', ecolor = 'navy', mec='navy',marker='o')

    uv_zB26, uv_sfrB26, uv_sfrB26errdn, uv_sfrB26errup = common.load_observation(obsdir, 'Global/SFRD_Compilation_UV_Bouwens26.dat', [0,1,2,3])
    uv_sfrB26errdn = uv_sfrB26 - uv_sfrB26errdn #make them relative errors
    uv_sfrB26errup = uv_sfrB26errup - uv_sfrB26

    uv_yobsB26 = uv_sfrB26 + np.log10(hobs/h0)

    #ax.errorbar(uv_zB26, uv_yobsB26, yerr=[uv_sfrB26errdn, uv_sfrB26errup], ls='None', mfc='None', ecolor = 'blue', mec='blue',marker='o', label = 'Bouwens+26 (UV)')

    ax.plot(zlist, np.log10(csfrd[:,0,0]), 'navy', linewidth=3, label ='UV (Shark v2.0 rr14)')
    ax.plot(zlist, np.log10(csfrd[:,0,1]), 'red', linewidth=3, label ='IR (Shark v2.0 rr14)')
    ax.plot(z_rr14steep, gcsfr_un, 'navy', linewidth=3, linestyle='dashed', label ='UV (Shark v2.0 rr14-steep)')
    ax.plot(z_rr14steep, gcsfr_at, 'red', linewidth=3, linestyle='dashed', label ='IR (Shark v2.0 rr14-steep)')

    print('#CSFRD global redshift CSFRD_unattenuated CSFRD_attenuated')
    for a,b,c,d,e,f,g in zip(zlist,np.log10(csfrd[:,0,0]), np.log10(csfrd[:,0,1]), np.log10(csfrd[:,1,0]), np.log10(csfrd[:,1,1]), np.log10(csfrd[:,0,0] - csfrd[:,1,0]), np.log10(csfrd[:,0,1] - csfrd[:,1,1])):
        print(a,b,c,d,e,f,g)
    common.prepare_legend(ax, ['navy','red','navy','red','darkgreen', 'darkred', 'blue'], loc=3) #bbox_to_anchor=(0.52, 0.47))

    ax = fig.add_subplot(132)
    plt.subplots_adjust(left=0.15)
    ax.text(1,-0.7,"Star formation in disks (L800)", fontsize=14)

    common.prepare_ax(ax, 0, 15, -6, -0.3, xtit, " ", locators=(1, 1, 1, 1))
    #plt.xscale('log')

    #ax.errorbar(zA23, yobsA23, yerr=[sfrA23errdn, sfrA23errup], ls='None', mfc='None', ecolor = 'darkgreen', mec='darkgreen',marker='s', label = 'Adams+24')
    #ax.errorbar(uv_zB26, uv_yobsB26, yerr=[uv_sfrB26errdn, uv_sfrB26errup], ls='None', mfc='None', ecolor = 'blue', mec='blue',marker='o', label = 'Bouwens+26 (UV)')

    ax.plot(zlist, np.log10(csfrd[:,1,0]), 'navy', linewidth=3, label ='UV (rr14)')
    ax.plot(zlist, np.log10(csfrd[:,1,1]), 'red', linewidth=3, label ='IR (rr14)')
    ax.plot(z_rr14steep, dcsfr_un, 'navy', linewidth=3, linestyle='dashed', label ='UV (rr14-steep)')
    ax.plot(z_rr14steep, dcsfr_at, 'red', linewidth=3, linestyle='dashed', label ='IR (rr14-steep)')

    #print('#CSFRD disks redshift CSFRD_unattenuated CSFRD_attenuated')
    #for a,b,c in zip(zlist,np.log10(csfrd[:,1,0]), np.log10(csfrd[:,1,1])):
    #    print(a,b,c)

    common.prepare_legend(ax, ['navy','red', 'navy','red','k','k','k'], loc=3) #bbox_to_anchor=(0.52, 0.47))

    ax = fig.add_subplot(133)
    plt.subplots_adjust(left=0.15)
    ax.text(3,-0.7,"Starbursts (L800)", fontsize=14)

    common.prepare_ax(ax, 0, 15, -6, -0.3, xtit, " ", locators=(1, 1, 1, 1))
    #plt.xscale('log')


    #ax.errorbar(zA23, yobsA23, yerr=[sfrA23errdn, sfrA23errup], ls='None', mfc='None', ecolor = 'darkgreen', mec='darkgreen',marker='s', label = 'Adams+24')
    #ax.errorbar(uv_zB26, uv_yobsB26, yerr=[uv_sfrB26errdn, uv_sfrB26errup], ls='None', mfc='None', ecolor = 'blue', mec='blue',marker='o', label = 'Bouwens+26 (UV)')

    ax.plot(zlist, np.log10(csfrd[:,0,0] - csfrd[:,1,0]), 'navy', linewidth=3, label ='UV (rr14)')
    ax.plot(zlist, np.log10(csfrd[:,0,1] - csfrd[:,1,1]), 'red', linewidth=3, label ='IR (rr14)')
    ax.plot(z_rr14steep, scsfr_un, 'navy', linewidth=3, linestyle='dashed', label ='UV (rr14-steep)')
    ax.plot(z_rr14steep, scsfr_at, 'red', linewidth=3, linestyle='dashed', label ='IR (rr14-steep)')

    #print('#CSFRD starbursts redshift CSFRD_unattenuated CSFRD_attenuated')
    #for a,b,c in zip(zlist,np.log10(csfrd[:,0,0] - csfrd[:,1,0]), np.log10(csfrd[:,0,1] - csfrd[:,1,1])):
    #    print(a,b,c)

    common.prepare_legend(ax, ['navy','red','navy','red', 'k','k','k'], loc=3) #bbox_to_anchor=(0.52, 0.47))
    plt.tight_layout()
    common.savefig(outdir, fig, "cosmic_sfr_from_lum_" + file_name + ".pdf")


def plot_mstar_mdust(plt, outdir, zlist, mdust_msm_redshift):

    #comparison with micro SURFS
    fig = plt.figure(figsize=(6,5.5))

    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(M_{\\rm dust}/M_{\\star})$"

    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15)

    common.prepare_ax(ax, 0, 7, -4, -1, xtit, ytit, locators=(1, 1, 1, 1))
    #plt.xscale('log')
  
    cols = ['Moccasin', 'PaleGoldenrod', 'Gold', 'DarkOrange', 'OrangeRed', 'LightCoral', 'Crimson', 'DarkRed', 'DarkViolet', 'Purple', 'Indigo']
    for j,b in enumerate(xmf):
        ind = np.where(mdust_msm_redshift[:,j,3] != 0)
        if(len(zlist[ind]) > 0):
           print(ind)
           x = zlist[ind]
           y = np.log10(mdust_msm_redshift[ind,j,3][0])
           ydn = np.log10(mdust_msm_redshift[ind,j,1][0])
           yup = np.log10(mdust_msm_redshift[ind,j,2][0])
   
           ax.fill_between(x,ydn, yup,facecolor=cols[j], alpha=0.4,interpolate=True) 
           ax.plot(x,y,linestyle='solid',color=cols[j], label = '%s'%str(b))
   
    common.prepare_legend(ax, cols, loc=1) #bbox_to_anchor=(0.52, 0.47))
    plt.tight_layout()
    common.savefig(outdir, fig, "Mdust_Mstar_redshift.pdf")



def prepare_data(hdf5_data, seds_lir, seds_bands, index, redshift, csfrd, csfrd_lums, mdust_msm_redshift):
   
    (h0, volh, sfr_disk, sfr_burst, mgasd, mgasb, mzd, mzb, msd, msb) = hdf5_data

    seds_total = seds_bands[1]
    lir_total = seds_lir[1]
    uv_total = seds_total[17,:] #band 17 is = FUV_Nathan a top-hat filter around 1500Angs 

    seds_disk = seds_bands[0]
    lir_disk = seds_lir[0]
    uv_disk = seds_disk[17,:] #band 17 is = FUV_Nathan a top-hat filter around 1500Angs 


    bin_it = functools.partial(us.wmedians, xbins=xmf)
    ind = np.where((uv_total < -1) & (uv_total > -35))
    sfr_total_uv = 10**((uv_total[ind] + 48.6) / (-2.5)) * ab_to_lsun / conv_uv_sfr #in Msun/yr

    ind = np.where(lir_total > 0)
    sfr_total_ir = lir_total[ind] * Lsun / conv_ir_sfr #in Msun/yr

    ind = np.where((uv_disk < -1) & (uv_disk > -35))
    sfr_disk_uv = 10**((uv_disk[ind] + 48.6) / (-2.5)) * ab_to_lsun / conv_uv_sfr #in Msun/yr

    ind = np.where(lir_disk > 0)
    sfr_disk_ir = lir_disk[ind] * Lsun / conv_ir_sfr #in Msun/yr

    (mdustd, DToM_MW) = dust_mass(mzd, mgasd, h0)
    (mdustb, DToM_MW) = dust_mass(mzb, mgasb, h0)

    ms_tot = (msd+ msb)/h0
    md_tot = mdustd + mdustb

    vol = volh/h0**3
    sfr_disk = sfr_disk/1e9/h0
    sfr_burst = sfr_burst/1e9/h0
    sfr_total = (sfr_disk + sfr_burst)
    ssfr = sfr_total/ms_tot

    csfrd[index,0,0] = np.sum(sfr_total)/vol
    csfrd[index,1,0] = np.sum(sfr_disk)/vol
    csfrd[index,2,0] = np.sum(sfr_burst)/vol


    csfrd_lums[index,0,0] = np.sum(sfr_total_uv) / vol
    csfrd_lums[index,0,1] = np.sum(sfr_total_ir) / vol
    csfrd_lums[index,1,0] = np.sum(sfr_disk_uv) / vol
    csfrd_lums[index,1,1] = np.sum(sfr_disk_ir) / vol

    #add perturbation
    sigma_sf = 0.3
    err = np.random.normal(0.0, sigma_sf, len(sfr_disk))
    csfrd[index,0,1] = np.sum(sfr_total*10**(err))/vol
    csfrd[index,1,1] = np.sum(sfr_disk*10**(err))/vol
    csfrd[index,2,1] = np.sum(sfr_burst*10**(err))/vol

    for j,b in enumerate(xmf):
        inb =  np.where((ms_tot >= 10**xmf_edges[j]) & (ms_tot < 10**xmf_edges[j+1]) & (ssfr >= 0)) # 0.2/(us.hubble_time(redshift) * 1e9))) #in yr^-1 ))
        if(len(ms_tot[inb]) > 99):
            mdms = md_tot[inb]/ms_tot[inb]
            mdust_msm_redshift[index,j,0] = np.mean(mdms)
            mdust_msm_redshift[index,j,3] = np.median(mdms)
            mdust_msm_redshift[index,j,1:3] = np.percentile(mdms, [16.0,84.0])

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    zlist= np.array([0, 0.1, 0.2029, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.76, 0.8, 0.9139, 1.0, 1.15, 1.3, 1.5, 1.65, 1.8, 2.0, 2.26, 2.51, 2.76, 3.0, 3.533, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 10.5, 11.0,11.5, 12, 12.5, 13.13, 13.5, 14.0, 14.5, 15.0]) #, 16.0, 17.0 ])

    file_name = "eagle-rr14"
    file_hdf5_sed = "Shark-SED-JWST-" + file_name + ".hdf5"

    plt = common.load_matplotlib()
    csfrd = np.zeros(shape = (len(zlist), 3, 2))
    csfrd_lums = np.zeros(shape = (len(zlist), 2, 2))

    mdust_msm_redshift   = np.zeros(shape = (len(zlist), len(xmf), 4))

    fields = {'galaxies': ('sfr_disk', 'sfr_burst', 'mgas_disk', 'mgas_bulge', 'mgas_metals_disk',
                           'mgas_metals_bulge', 'mstars_disk', 'mstars_bulge')}

    fields_sed = {'SED/lir_dust': ('disk','total'),}
    fields_seds_bands = {'SED/ab_dust':('disk','total'),}

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        seds_lir = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed, subvols, file_hdf5_sed)
        seds_bands = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_seds_bands, subvols, file_hdf5_sed)

        mass = prepare_data(hdf5_data, seds_lir, seds_bands, index, zlist[index], csfrd, csfrd_lums, mdust_msm_redshift)

        h0 = hdf5_data[0]
        volh = hdf5_data[1]

    plot_SFRD(plt, outdir, obsdir, csfrd, zlist, h0)
    plot_SFRD_from_lum(plt, outdir, obsdir, csfrd_lums, zlist, h0, file_name)

    #plot_mstar_mdust(plt, outdir, zlist, mdust_msm_redshift)
    #np.savetxt("MstarMdust_Redshift_mean.txt", mdust_msm_redshift[:,:,0])
    #np.savetxt("MstarMdust_Redshift_16th.txt", mdust_msm_redshift[:,:,1])
    #np.savetxt("MstarMdust_Redshift_84th.txt", mdust_msm_redshift[:,:,2])
    #np.savetxt("MstarMdust_Redshift_median.txt", mdust_msm_redshift[:,:,3])


if __name__ == '__main__':
    main(*common.parse_args())
