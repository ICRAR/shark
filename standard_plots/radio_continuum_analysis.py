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

mslow = 7.0
msupp = 12.5
dms = 0.25
msbins = np.arange(mslow, msupp, dms)
xmfs = msbins + dms/2.0


mbhlow = 6.0
mbhupp = 10.
dmbh = 1
mbhbins = np.arange(mbhlow, mbhupp, dmbh)
xmbhf = mbhbins + dmbh/2.0


zsun = 0.0189
MpctoKpc = 1e3
PI       = 3.141592654
h        = 6.6261e-27 #cm2 g s-1
d10pc = 3.086e+19 #10 parsecs in cm
dfac = 4 * PI * d10pc**2

#choose radio continuum model
radio_model = "Bressan02"
c_speed = 299792458.0 #in m/s
c_speed_cm = c_speed * 1e2 #in cm/s
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


def ionising_photons(m, wave):
    #m is a vector of AB absolute magnitudes in a band with central wavelength wave
    #wavelength input has to be in angstrom
    
    wave_m = wave * 1e-10 #wavelength in m
    wave_cm = wave_m  * 1e2 #wavelength in cm
    freq = c_speed / wave_m #Hz
    hc =  h * (c_speed * 1e2) #h*c in cgs

    lum = 10.0**((m + 48.6) / (-2.5)) * dfac * freq * wave_cm #we want to convert from Luminosity to int(lambda*Lum_lambda*dlambda)
    Q = lum / hc #rate of ionising photons in s^-1.
  
    return Q

def freefree_lum(Q, nu):

    #Q is the rate of ionising photos in s^-1
    #nu is the frequency in GHz
    #output in erg/s/Hz

    T = 1e4 #temperature in K
    lum = Q/6.3e32 * (T/1e4)**0.45 * (nu)**(-0.1)


    return lum

def synchrotron_lum(SFR, nu):
   
    #SFR in Msun/yr
    #nu is the frequency in GHz
    #output in erg/s/Hz

    ENT = 1.44 
    ESNR = 0.06 * ENT
    alpha = -0.8

    T = 1e4
    EM = 6.5 #pc * cm**-6; this comes from the observations of H. V. Cane (1979-11) Spectra of the non-thermal radio radiation from the galactic polar regions. MNRAS 189, pp. 465–478. Cited by: §4.3.2, who found that τ≈1 at ν≈2 MHz.
    tau = (T/1e4)**(-1.35) * (nu / 1.4)**(-2.1) * EM / 6e6
    comp1 = ESNR * (nu / 1.49)**(-0.5) + ENT * (nu / 1.49)**(alpha) * np.exp(-tau)
    nuSNCC = SFR * 0.011148
    lum = comp1 * 1e30 * nuSNCC

    return lum


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

def prepare_data(hdf5_data, seds_nod, seds, lir, seds_lir_bc, index, model_dir, snapshot, filters, hist_sf_dale14, hist_sf_obi17, hist_agn, hist_agn_bh, redshift):


    bin_it = functools.partial(us.wmedians, xbins=xmf)
    total_mags_nod = seds_nod[4]
    total_mags = seds[4]

    lir_bc_cont = seds_lir_bc[2]
    Lum_radio_Viperfish = 10**((total_mags[9:16,:] + 48.6)/(-2.5)) * dfac #erg/s/Hz

    ion_mag = total_mags_nod[1,:]
    q_ionis = ionising_photons(ion_mag, 912.0) #in s^-1

    # Unpack data
    (h0, volh, mdisk, mbulge, sfrd, sfrb, idgal, mbh, macc_hh, macc_sb, mgd, mgb, typeg) = hdf5_data

    temp_bc = 57.0
    temp_diff = 22.0
    temp_eff =  temp_bc * lir_bc_cont + temp_diff * (1.0 - lir_bc_cont)


    mbh = mbh/h0
    macc_bh = (macc_hh + macc_sb)/h0/1e9 #in Msun/yr
    (Ljet_ADAF, Ljet_td, mdot_norm) = radio_luminosity_agn(mbh, macc_bh)

    h0log = np.log10(float(h0))

    vol = volh/h0**3

    sfr = sfrd + sfrb

    #select galaxies with Mstar > 0
    ind = np.where((mdisk + mbulge  > 0) | (sfr > 0))
    ms = (mdisk[ind] + mbulge[ind])/h0 #in Msun
    sfr = sfr[ind]/h0/1e9 #in Msun/yr
    typeg = typeg[ind]

    #calculate main sequence
    ind = np.where((ms > 3e8) & (ms <1e10) & (typeg ==0) & (sfr > 0))
    ms_fit = np.polyfit(np.log10(ms[ind]), np.log10(sfr[ind]), 1)
    delta_ms = np.log10(sfr) - (ms_fit[0] * np.log10(ms) + ms_fit[1])
    print("Delta main seq", delta_ms, np.median(delta_ms))

 
    selection_freq = (8.4, 5.0, 3.0, 1.4, 0.61, 0.325, 0.15) #GHz

    lum_radio = np.zeros(shape = (len(selection_freq), len(q_ionis)))
    lum_ratio = np.zeros(shape = (len(selection_freq), len(q_ionis)))
    lum_radio_agn = np.zeros(shape = (len(selection_freq), len(mbh)))

    for i, nu in enumerate(selection_freq):
        lum_radio[i,:] = freefree_lum(q_ionis[:], nu) + synchrotron_lum(sfr[:], nu)
        lum_ratio[i,:] = freefree_lum(q_ionis[:], nu) / lum_radio[i,:]
        lum_radio_agn[i,:] = radio_luminosity_per_freq(Ljet_ADAF[:], Ljet_td[:], mdot_norm[:], mbh[:], nu)

    max_ff = max(lum_ratio[3,:] * lum_radio[3,:])
    sfr_derived = lum_ratio[3,:] * lum_radio[3,:] / (2.17e27) * (1.4)**0.1
    max_sfr_derived = max_ff / (2.17e27) * (1.4)**0.1

    lir_total = lir[1] #total dust luminosity
    qIR_dale14 = np.log10(lir_total[0,:]*Lsunwatts/3.75e12) - np.log10(Lum_radio_Viperfish[3,:]/1e7)
    qIR_bressan = np.log10(lir_total[0,:]*Lsunwatts/3.75e12) - np.log10(lum_radio[3,:]/1e7)
 
    ind = np.where(lir_total[0,:] > 1e7)
    print("Median qIR for Dale14 and Bressan", np.median(qIR_dale14[ind]), np.median(qIR_bressan[ind])," at redshift", redshift)

    ind= np.where((ms > 1e8) & (ms < 1e9) & (temp_eff[0,:] > 50) & (typeg ==0))
    print("Median sfr, ms, qIR", np.median(sfr[ind]), np.median(ms[ind]), np.median(qIR_bressan[ind]))

    ind= np.where((ms > 1e8) & (ms < 1e9) & (qIR_bressan  > 2.5 ) & (typeg == 0))
    print("Median sfr, ms, qIR", np.median(sfr[ind]), np.median(ms[ind]), np.median(qIR_bressan[ind]))

    ind = np.where(sfr > 300)
    print(lir_total[0,ind])

    ran_err = np.random.normal(0.0, 0.4, len(sfr))

    ind = np.where(Lum_radio_Viperfish[3,:]/1e7 > 1e17)
    H, _ = np.histogram(np.log10(Lum_radio_Viperfish[3,ind]/1e7),bins=np.append(mbins,mupp))
    hist_sf_dale14[index,0,:] = hist_sf_dale14[index,0,:] + H
    H, _ = np.histogram(np.log10(Lum_radio_Viperfish[3,ind]/1e7) + ran_err[ind],bins=np.append(mbins,mupp))
    hist_sf_dale14[index,1,:] = hist_sf_dale14[index,1,:] + H

    ind = np.where(lum_radio[3,:]/1e7 > 1e17)
    H, _ = np.histogram(np.log10(lum_radio[3,ind]/1e7),bins=np.append(mbins,mupp))
    hist_sf_obi17[index,0,:] = hist_sf_obi17[index,0,:] + H
    H, _ = np.histogram(np.log10(lum_radio[3,ind]/1e7) + ran_err[ind],bins=np.append(mbins,mupp))
    hist_sf_obi17[index,1,:] = hist_sf_obi17[index,1,:] + H

    ind = np.where(lum_radio_agn[3,:]/1e7 > 1e17)
    H, _ = np.histogram(np.log10(lum_radio_agn[3,ind]/1e7),bins=np.append(mbins,mupp))
    hist_agn[index,:] = hist_agn[index,:] + H

    for j,m in enumerate(xmbhf):
            ind = np.where((lum_radio_agn[3,:]/1e7 > 1e17) & (np.log10(mbh) >= m - dmbh/2.0) & (np.log10(mbh) < m + dmbh/2.0))
            H, _ = np.histogram(np.log10(lum_radio_agn[3,ind]/1e7),bins=np.append(mbins,mupp))
            hist_agn_bh[index,j,:] = hist_agn_bh[index,j,:] + H


    hist_sf_dale14[index,:] = hist_sf_dale14[index,:]/vol/dm
    hist_sf_obi17[index,:] = hist_sf_obi17[index,:]/vol/dm
    hist_agn[index,:] = hist_agn[index,:]/vol/dm
    hist_agn_bh[index,:] = hist_agn_bh[index,:]/vol/dm

    #dividing by 1e7 the luminosities to output them in W/Hz 
    return (lum_radio/1e7, Lum_radio_Viperfish/1e7, lum_ratio, lum_radio_agn/1e7, ms, sfr, vol, h0, qIR_bressan, temp_eff, typeg, delta_ms)

def plot_comparison_radio_lums(plt, outdir, obsdir, LBressan, LViperfish, Lratio, ms, sfr, qIR_bressan, Tdust, typeg, delta_ms, filters, redshift):

    bin_it = functools.partial(us.wmedians, xbins=xmf)
    bin_it_sm = functools.partial(us.wmedians, xbins=xmfs)

    Tdust = Tdust[0]

    subplots = (331, 332, 333, 334, 335, 336, 337)
    labels = filters

    fig = plt.figure(figsize=(10,10))
    ytit = "$\\rm log_{10} (L_{\\rm Bressan02}/W Hz^{-1})$"
    xtit = "$\\rm log_{10} (L_{\\rm Dale14}/W Hz^{-1})$"
    xmin, xmax, ymin, ymax = 19, 25, 19, 25
    xleg = xmin + 0.15 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    for i, subp in enumerate(subplots):
        ax = fig.add_subplot(subp)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))
        ax.text(xleg, yleg, labels[i])
 
        ind = np.where((ms > 1e8) & (sfr > 0) & (LViperfish[i,:] > 1e19))
        x = np.log10(LViperfish[i,ind])
        y = np.log10(LBressan[i,ind])
        meds =  bin_it(x=x[0], y=y[0])
        z = np.log10(Lratio[i,ind])
        im = ax.hexbin(x[0], y[0], z[0], xscale='linear', yscale='linear', gridsize=(20,20), cmap='magma')
        #ax.plot(x[0], y[0], 'ko')
        #cbar_ax = fig.add_axes([0.86, 0.15, 0.025, 0.7])
        cbar = fig.colorbar(im)#, cax=cbar_ax)
        cbar.ax.set_ylabel('$\\rm log_{10}(L_{\\rm ff}/L_{\\rm syn}(Bressan)$)')
        x = [19, 25]
        ax.plot(x, x, linestyle='dotted',color='yellow')

        ind = np.where(meds[0,:] != 0)
        x = xmf[ind]
        y = meds[0,ind] 
        yerrdn = meds[1,ind]
        yerrup = meds[2,ind]
        ax.errorbar(x, y[0], yerr=[yerrdn[0], yerrup[0]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')
    #common.prepare_legend(ax, cols, loc=3)

    plt.tight_layout()

    common.savefig(outdir, fig, 'radio_lum_comparison_models_z'+redshift+'.pdf')


    subplots = (211, 212)
    labels = ['Dale14','Bressan02']

    fig = plt.figure(figsize=(4.5,8))
    ytit = "$\\rm log_{10} (SFR/M_{\\odot}\\, yr^{-1})$"
    xtit = "$\\rm log_{10} (L_{\\rm 1.4GHz}/W Hz^{-1})$"
    xmin, xmax, ymin, ymax = 19, 25, -3, 3
    xleg = xmin + 0.15 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    for i, subp in enumerate(subplots):
        if i == 0:
           Lr = LViperfish[3,:]
        else:
           Lr = LBressan[3,:]
        ax = fig.add_subplot(subp)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))
        ax.text(xleg, yleg, labels[i])
 
        ind = np.where((ms > 1e8) & (sfr > 0) & (Lr[:] > 1e19))
        x = np.log10(Lr[ind])
        y = np.log10(sfr[ind])
        meds =  bin_it(x=x, y=y)
        im = ax.hexbin(x, y, xscale='linear', yscale='linear', gridsize=(20,20), cmap='magma', mincnt=4)

        ind = np.where(meds[0,:] != 0)
        x = xmf[ind]
        y = meds[0,ind] 
        yerrdn = meds[1,ind]
        yerrup = meds[2,ind]
        ax.errorbar(x, y[0], yerr=[yerrdn[0], yerrup[0]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')
    #common.prepare_legend(ax, cols, loc=3)

    plt.tight_layout()

    common.savefig(outdir, fig, 'radio_lum_sfr_models_z'+redshift+'.pdf')


    fig = plt.figure(figsize=(4.5,8))
    xtit = "$\\rm log_{10} (M_{\\star}/M_{\\odot})$"
    ytit = "$\\rm q_{\\rm IR}$"
    xmin, xmax, ymin, ymax = 8, 12, 0, 3
    xleg = xmin + 0.15 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))
    ax.text(xleg, yleg, 'Bressan model')

    ind = np.where((qIR_bressan > -1) & (ms > 1e8) & (qIR_bressan < 5) & (abs(delta_ms) < 0.3))
    x = np.log10(ms[ind])
    y = qIR_bressan[ind]
    z = Tdust[ind]
    meds =  bin_it_sm(x=x, y=y)
    im = ax.hexbin(x, y, xscale='linear', yscale='linear', gridsize=(20,20), cmap='magma', mincnt=4)

    ind = np.where(meds[0,:] != 0)
    xm = xmfs[ind]
    ym = meds[0,ind]
    yerrdn = meds[1,ind]
    yerrup = meds[2,ind]
    ax.errorbar(xm, ym[0], yerr=[yerrdn[0], yerrup[0]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')
    #common.prepare_legend(ax, cols, loc=3)

    ax = fig.add_subplot(212)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))

    im = ax.hexbin(x, y, z, xscale='linear', yscale='linear', gridsize=(20,20), cmap='magma', mincnt=4)
    cbar = fig.colorbar(im)#, cax=cbar_ax)
    cbar.ax.set_ylabel('$\\rm T_{\\rm dust}/K$')

    ax.errorbar(xm, ym[0], yerr=[yerrdn[0], yerrup[0]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')
    #common.prepare_legend(ax, cols, loc=3)


    plt.tight_layout()

    common.savefig(outdir, fig, 'qIR_sfr_z'+redshift+'.pdf')

def plot_radio_lf(plt, output_dir, obs_dir, hist_sf_dale14, hist_sf_obi17, hist_agn, hist_agn_bh, h0):

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
        xmin, xmax, ymin, ymax = 19, 26, -7, -1
        xleg = xmax - 0.15 * (xmax - xmin)
        yleg = ymax - 0.1 * (ymax - ymin)
        ax = fig.add_subplot(111)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))
        ax.text(xleg, yleg, zlabel[i]) 
        ind = np.where(hist_sf_obi17[i,0,:] != 0)
        xp = xmf[ind]
        yp = np.log10(hist_sf_obi17[i,0,ind])
        ax.plot(xp, yp[0], linestyle='solid', color='red', label='Shark + Obi17')
        ind = np.where(hist_sf_obi17[i,1,:] != 0)
        xp = xmf[ind]
        yp = np.log10(hist_sf_obi17[i,1,ind])
        ax.plot(xp, yp[0], linestyle='dotted', color='red', label='+0.3dex err')

        ind = np.where(hist_sf_dale14[i,0,:] != 0)
        xp = xmf[ind]
        yp = np.log10(hist_sf_dale14[i,0,ind])
        ax.plot(xp, yp[0], linestyle='solid', color='blue', label='Shark + Dale14')
        ind = np.where(hist_sf_dale14[i,1,:] != 0)
        xp = xmf[ind]
        yp = np.log10(hist_sf_dale14[i,1,ind])
        ax.plot(xp, yp[0], linestyle='dotted', color='blue', label='+0.3dex err')

   
        rad_obs(ax, zf, 'SF') 
     
        common.prepare_legend(ax, ['red', 'red','blue', 'blue','grey'], loc=3)
        plt.tight_layout()
        common.savefig(output_dir, fig, 'radio_lum_function_' + zf + '.pdf')
    
    
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
    xmin, xmax, ymin, ymax = 19, 26, -7, -1
    xleg = xmax - 0.15 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)
    subp = (231, 232, 233, 234, 235, 236)

    for i, zf in enumerate(zfiles):
        ax = fig.add_subplot(subp[i])
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))
        ax.text(xleg, yleg, zlabel[i]) 
        ind = np.where(hist_sf_obi17[i,0,:] != 0)
        xp = xmf[ind]
        yp = np.log10(hist_sf_obi17[i,0,ind])
        ax.plot(xp, yp[0], linestyle='solid', color='red', label='Shark + Obi17')
        ind = np.where(hist_sf_obi17[i,1,:] != 0)
        xp = xmf[ind]
        yp = np.log10(hist_sf_obi17[i,1,ind])
        ax.plot(xp, yp[0], linestyle='dotted', color='red', label='+0.3dex err')

        ind = np.where(hist_sf_dale14[i,0,:] != 0)
        xp = xmf[ind]
        yp = np.log10(hist_sf_dale14[i,0,ind])
        ax.plot(xp, yp[0], linestyle='solid', color='blue', label='Shark + Dale14')
        ind = np.where(hist_sf_dale14[i,1,:] != 0)
        xp = xmf[ind]
        yp = np.log10(hist_sf_dale14[i,1,ind])
        ax.plot(xp, yp[0], linestyle='dotted', color='blue', label='+0.3 dex')
  
        rad_obs(ax, zf, 'SF') 
     
        common.prepare_legend(ax, ['red', 'red', 'blue', 'blue', 'grey'], loc=3)
    plt.tight_layout()
    common.savefig(output_dir, fig, 'radio_lum_function_allz.pdf')
    
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
        ax.plot(xp, yp[0], linestyle='solid', color='black', label='Shark (all AGN)')

        rad_obs(ax, zf, 'AGN') 
     
        common.prepare_legend(ax, ['black', 'grey'], loc=3)
    plt.tight_layout()
    common.savefig(output_dir, fig, 'radio_lum_function_AGN_allz_v2.pdf')
 
def main(model_dir, output_dir, redshift_table, subvols, obs_dir):

    #zlist = np.arange(2,10,0.25)
    #zlist = (0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 0, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 9, 10)
    filters = ('8.4GHz', '5GHz', '3GHz', '1.4GHz', '610MHz', '325MHz', '150MHz') 

    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5" #"Shark-SED-eagle-rr14-radio-only.hdf5"
    #(0): "z_SDSS", "Band_ionising_photons", "FUV_Nathan", "Band9_ALMA",
    #(4): "Band8_ALMA", "Band7_ALMA", "Band6_ALMA", "Band4_ALMA", "Band3_ALMA",
    #(9): "BandX_VLA", "BandC_VLA", "BandS_VLA", "BandL_VLA", "Band_610MHz",
    #(14): "Band_325MHz", "Band_150MHz"

    #199 188 159 131 113 100 88 79 70 63 57 51
    zlist = [0, 1.0, 2.00391410007239, 3.0191633709527, 3.95972701662501, 5.02220991014863, 5.96592270612165, 7.05756323172746, 8.0235605165086, 8.94312532315157, 9.95650268434316] #9.95655] #0.194739, 0.254144, 0.359789, 0.450678, 0.8, 0.849027, 0.9, 1.20911, 1.28174, 1.39519, 1.59696, 2.00392, 2.47464723643932, 2.76734390952347, 3.01916, 3.21899984389701, 3.50099697082904, 3.7248038025221, 3.95972, 4.465197621546, 4.73693842543988] #[5.02220991014863, 5.52950356184419, 5.96593, 6.55269895697227, 7.05756323172746, 7.45816170313544, 8.02352, 8.94312532315157, 9.95655]
    #[0.016306640039433, 0.066839636933135, 0.084236502339783, 0.119886040396529, 0.138147164704691, 0.175568857770275, 0.214221447279112, 0.23402097095238, 0.274594901875312, 0.316503156974571]

    znames = ['0', '0p2', '0p9', '2', '3', '4', '5', '6', '7', '8', '9', '10']
    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge','sfr_disk','sfr_burst','id_galaxy',
                           'm_bh', 'bh_accretion_rate_hh', 'bh_accretion_rate_sb', 'mgas_disk', 'mgas_bulge', 'type')}
   
    fields_sed_nod = {'SED/ab_nodust': ('bulge_d','bulge_m','bulge_t','disk','total')}
    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total')}
    fields_lir = {'SED/lir_dust': ('disk','total')}
    fields_seds_bc = {'SED/lir_dust_contribution_bc': ('disk','bulge_t','total'),}

    hist_sf_dale14  = np.zeros(shape = (len(zlist), 2, len(mbins)))
    hist_sf_obi17  = np.zeros(shape = (len(zlist), 2, len(mbins)))
    hist_agn  = np.zeros(shape = (len(zlist), len(mbins)))

    hist_agn_bh  = np.zeros(shape = (len(zlist), len(mbhbins),len(mbins)))


    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        seds_nod = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_nod, subvols, file_hdf5_sed)
        seds = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed)
        lir = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_lir, subvols, file_hdf5_sed)
        seds_lir_bc = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_seds_bc, subvols, file_hdf5_sed)

        (LBressan, LViperfish, Lratio, LAGN, ms, sfr, vol, h0, qIR_bressan, tdust, typeg, delta_ms) = prepare_data(hdf5_data, seds_nod, seds, lir, seds_lir_bc, index, model_dir, snapshot, filters, hist_sf_dale14, hist_sf_obi17, hist_agn, hist_agn_bh, zlist[index])
        plot_comparison_radio_lums(plt, output_dir, obs_dir, LBressan, LViperfish, Lratio, ms, sfr, qIR_bressan, tdust, typeg, delta_ms, filters, znames[index])

    plot_radio_lf(plt, output_dir, obs_dir, hist_sf_dale14, hist_sf_obi17, hist_agn, hist_agn_bh, h0)

if __name__ == '__main__':
    main(*common.parse_args())
