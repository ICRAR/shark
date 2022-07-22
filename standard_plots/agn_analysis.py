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
import os
import h5py

import common
import utilities_statistics as us

##################################
# Constants
mlow = 18.0
mupp = 28.0
dm = 0.25
mbins = np.arange(mlow, mupp, dm)
xmf = mbins + dm/2.0

mbhlow = 6.0
mbhupp = 10.0
dmbh = 1
mbhbins = np.arange(mbhlow, mbhupp, dmbh)
xmbhf = mbhbins + dmbh/2.0


zsun = 0.0189
MpctoKpc = 1e3
PI       = 3.141592654
h        = 6.6261e-27 #cm2 g s-1
d10pc = 3.086e+19 #10 parsecs in cm
dfac = 4 * PI * d10pc**2

#model of Mattson et al. (2014) for the dependence of the dust-to-metal mass ratio and metallicity X/H.
corrfactor_dm = 2.0
polyfit_dm = [ 0.00544948, 0.00356938, -0.07893235,  0.05204814,  0.49353238] 

#choose radio continuum model
radio_model = "Bressan02"
c_speed = 299792458.0 #in m/s
c_speed_cm = 299792458.0 * 1e2 #in cm/s
Lsunwatts = 3.846e26

#define parameters associated to AGN luminosity
alpha_adaf = 0.1
delta = 0.0005
beta = 1 - alpha_adaf / 0.55
eta_edd = 4
mcrit_nu = 0.001 * (delta / 0.0005) * (1 - beta) / beta * alpha_adaf**2
spin = 0.1
A_adaf = 4e-4
A_td = A_adaf/100



def radio_luminosity_agn(mbh, macc):

    #input mbh has to be in Msun
    #input macc has to be in Msun/yr
    #output luminosity in erg/s
   
    Ljet_ADAF = np.zeros ( shape = len(mbh))
    Ljet_td = np.zeros ( shape = len(mbh))
    mdot_norm = np.zeros ( shape = len(mbh))
    Ledd = np.zeros ( shape = len(mbh))

    ind = np.where(mbh > 0)
    Ledd[ind] = 1.28e46 * (mbh[ind]/1e8) #in erg/s

    ind = np.where((mbh > 0) & (macc > 0))
    mdot_edd = Ledd[ind] / (0.1 * c_speed_cm**2) * 1.586606334841629e-26 #in Msun/yr
    mdot_norm[ind] = macc[ind] / mdot_edd

    Ljet_ADAF[ind] = 2e45 * (mbh[ind]/1e9) * (mdot_norm[ind]/0.01) * spin**2 #in erg/s
    Ljet_td[ind] = 2.5e43 * (mbh[ind]/1e9)**1.1 * (mdot_norm[ind]/0.01)**1.2 * spin**2 #in erg/s

    ind = np.where(mdot_norm > eta_edd)
    mdot_norm[ind] = eta_edd

    return (Ljet_ADAF, Ljet_td, mdot_norm)

def radio_luminosity_per_freq(Ljet_ADAF, Ljet_td, mdot_norm, mbh, nu):

    #here nu has to be rest-frame in GHz
    lum_1p4GHz_adaf = A_adaf * (mbh / 1e9 * mdot_norm / 0.01)**0.42 * Ljet_ADAF
    lum_1p4GHz_td = A_td * (mbh / 1e9)**0.32 * (mdot_norm / 0.01)**(-1.2) * Ljet_td

    freq_hz = nu * 1e9
    lum_nu = (lum_1p4GHz_adaf + lum_1p4GHz_td) * (nu / 1.4)**(-0.7) / freq_hz #in erg/s/Hz

    return lum_nu

def prepare_data(hdf5_data, index, model_dir, snapshot, filters, hist_agn, hist_agn_bh):


    bin_it = functools.partial(us.wmedians, xbins=xmf)

    # Unpack data
    (h0, volh, mdisk, mbulge, sfrd, sfrb, idgal, mbh, macc_hh, macc_sb) = hdf5_data

    mbh = mbh/h0
    macc_bh = (macc_hh + macc_sb)/h0/1e9 #in Msun/yr
    (Ljet_ADAF, Ljet_td, mdot_norm) = radio_luminosity_agn(mbh, macc_bh)

    h0log = np.log10(float(h0))

    vol = volh/h0**3

    sfr = sfrd + sfrb

    #select galaxies with Mstar > 0
    ind = np.where(mdisk + mbulge > 0)
    ms = (mdisk[ind] + mbulge[ind])/h0 #in Msun
    sfr = sfr[ind]/h0/1e9 #in Msun/yr
    
    selection_freq = (8.4, 5.0, 3.0, 1.4, 0.61, 0.325, 0.15) #GHz

    lum_radio_agn = np.zeros(shape = (len(selection_freq), len(mbh)))

    for i, nu in enumerate(selection_freq):
        lum_radio_agn[i,:] = radio_luminosity_per_freq(Ljet_ADAF[:], Ljet_td[:], mdot_norm[:], mbh[:], nu)

    ind = np.where(lum_radio_agn[3,:]/1e7 > 1e17)
    H, _ = np.histogram(np.log10(lum_radio_agn[3,ind]/1e7),bins=np.append(mbins,mupp))
    hist_agn[index,:] = hist_agn[index,:] + H

    for j,m in enumerate(xmbhf):
            ind = np.where((lum_radio_agn[3,:]/1e7 > 1e17) & (np.log10(mbh) >= m - dmbh/2.0) & (np.log10(mbh) < m + dmbh/2.0))
            H, _ = np.histogram(np.log10(lum_radio_agn[3,ind]/1e7),bins=np.append(mbins,mupp))
            hist_agn_bh[index,j,:] = hist_agn_bh[index,j,:] + H

    hist_agn[index,:] = hist_agn[index,:]/vol/dm
    hist_agn_bh[index,:] = hist_agn_bh[index,:]/vol/dm

    return (lum_radio_agn/1e7, ms, sfr, vol, h0)

def plot_radio_lf(plt, output_dir, obs_dir, hist_agn, hist_agn_bh, h0):

    def rad_obs(ax, z, Galaxies):
        lm, p, dpdn, dpup = common.load_observation(obs_dir, 'lf/lf1p4GHz_' + Galaxies+'_' + z +'_Bonato21.dat', [0,1,2,3])
        hobs = 0.7
        xobs = lm + 2.0 * np.log10(hobs/h0)
        yobs = p - 3.0 * np.log10(hobs/h0) 
        ax.errorbar(xobs, yobs, yerr=[dpdn, dpup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Bonato et al. (2020)")

    zfiles = ['z0', 'z1', 'z2', 'z3', 'z4', 'z5']
    zlabel = ['z=0', 'z=1', 'z=2', 'z=3', 'z=4', 'z=5']
    colbh = ['DarkBlue','DarkCyan','Gold','OrangeRed','DarkRed']
    labelbh = ['[6-7]','[7,8]','[8-9]','[9-10]', '>10']

    for i, zf in enumerate(zfiles):
   
        fig = plt.figure(figsize=(5,5))
        ytit = "$\\rm log_{10} (\\phi/Mpc^{-3} dex^{-1})$"
        xtit = "$\\rm log_{10} (L_{\\rm 1.4GHz}/W Hz^{-1})$"
        xmin, xmax, ymin, ymax = 19, 28.5, -7, -1
        xleg = xmax - 0.15 * (xmax - xmin)
        yleg = ymax - 0.1 * (ymax - ymin)
        ax = fig.add_subplot(111)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))
        ax.text(xleg, yleg, zlabel[i]) 
   
        ind = np.where(hist_agn[i,:] != 0)
        xp = xmf[ind]
        yp = np.log10(hist_agn[i,ind])
        ax.plot(xp, yp[0], linestyle='solid', color='black', label='Shark + Griffin19')
        for j,m  in enumerate(xmbhf):
            ind = np.where(hist_agn_bh[i,j,:] != 0)
            xp = xmf[ind]
            yp = np.log10(hist_agn_bh[i,j,ind])
            ax.plot(xp, yp[0], linestyle='solid', color=colbh[j], label=labelbh[j])

        rad_obs(ax, zf, 'AGN')    

        common.prepare_legend(ax, ['black'], loc=3)
        plt.tight_layout()
        common.savefig(output_dir, fig, 'radio_lum_function_AGN_' + zf + '.pdf')


    fig = plt.figure(figsize=(12,8))
    ytit = "$\\rm log_{10} (\\phi/Mpc^{-3} dex^{-1})$"
    xtit = "$\\rm log_{10} (L_{\\rm 1.4GHz}/W Hz^{-1})$"
    xmin, xmax, ymin, ymax = 19, 28.5, -8, -2
    xleg = xmax - 0.15 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)
    subp = (231, 232, 233, 234, 235, 236)

    for i, zf in enumerate(zfiles):
        ax = fig.add_subplot(subp[i])
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))
        ax.text(xleg, yleg, zlabel[i]) 
        ind = np.where(hist_agn[i,:] != 0)
        xp = xmf[ind]
        yp = np.log10(hist_agn[i,ind])
        ax.plot(xp, yp[0], linestyle='solid', color='black', label='Shark')
        for j,m  in enumerate(xmbhf):
            ind = np.where(hist_agn_bh[i,j,:] != 0)
            xp = xmf[ind]
            yp = np.log10(hist_agn_bh[i,j,ind])
            ax.plot(xp, yp[0], linestyle='solid', color=colbh[j], label=labelbh[j])

        rad_obs(ax, zf, 'AGN') 
     
        common.prepare_legend(ax, ['black', 'grey'], loc=3)
    plt.tight_layout()
    common.savefig(output_dir, fig, 'radio_lum_function_AGN_allz.pdf')
 

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):

    #zlist = np.arange(2,10,0.25)
    #zlist = (0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 0, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 9, 10)
    filters = ('8.4GHz', '5GHz', '3GHz', '1.4GHz', '610MHz', '325MHz', '150MHz') 

    #199 188 159 131 113 100 88 79 70 63 57 51
    zlist = [0, 1.0, 2.00391410007239, 3.0191633709527, 3.95972701662501, 5.02220991014863, 5.96592270612165, 7.05756323172746, 8.0235605165086, 8.94312532315157, 9.95650268434316] #, 0.194738848008908, 0.909822023685613, 2.00391410007239, 3.0191633709527, 3.95972701662501, 5.02220991014863, 5.96592270612165, 7.05756323172746, 8.0235605165086, 8.94312532315157, 9.95650268434316] #9.95655] #0.194739, 0.254144, 0.359789, 0.450678, 0.8, 0.849027, 0.9, 1.20911, 1.28174, 1.39519, 1.59696, 2.00392, 2.47464723643932, 2.76734390952347, 3.01916, 3.21899984389701, 3.50099697082904, 3.7248038025221, 3.95972, 4.465197621546, 4.73693842543988] #[5.02220991014863, 5.52950356184419, 5.96593, 6.55269895697227, 7.05756323172746, 7.45816170313544, 8.02352, 8.94312532315157, 9.95655]
    #[0.016306640039433, 0.066839636933135, 0.084236502339783, 0.119886040396529, 0.138147164704691, 0.175568857770275, 0.214221447279112, 0.23402097095238, 0.274594901875312, 0.316503156974571]

    znames = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge','sfr_disk','sfr_burst','id_galaxy',
                           'm_bh', 'bh_accretion_rate_hh', 'bh_accretion_rate_sb')}
    hist_agn  = np.zeros(shape = (len(zlist), len(mbins)))

    hist_agn_bh  = np.zeros(shape = (len(zlist), len(mbhbins),len(mbins)))

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        (LAGN, ms, sfr, vol, h0) = prepare_data(hdf5_data, index, model_dir, snapshot, filters, hist_agn, hist_agn_bh)
           
    plot_radio_lf(plt, output_dir, obs_dir, hist_agn, hist_agn_bh, h0)

if __name__ == '__main__':
    main(*common.parse_args())
