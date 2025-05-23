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
"""HMF plots"""

import numpy as np
import os

import common


##################################

# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3
c_speed = 299792458.0 #m/s

#choose radio continuum model
c_speed_cm = c_speed * 1e2
Lsunwatts = 3.846e26


MpctoKpc = 1e3
PI       = 3.141592654
h        = 6.6261e-27 #cm2 g s-1
d10pc = 3.086e+19 #10 parsecs in cm
dfac = 4 * PI * d10pc**2

# Mass function initialization

mlow = 35
mupp = 45
dm = 0.25
mbins = np.arange(mlow,mupp,dm)
xlf   = mbins + dm/2.0

def ionising_photons(m, wave):
    #m is a vector of AB absolute magnitudes in a band with central wavelength wave
    #wavelength input has to be in angstrom

    Q = np.zeros(shape = len(m))
    
    wave_m = wave * 1e-10 #wavelength in m
    wave_cm = wave_m  * 1e2 #wavelength in cm
    freq = c_speed / wave_m #Hz
    hc =  h * (c_speed * 1e2) #h*c in cgs

    ind = np.where(m != -999)
    lum = 10.0**((m[ind] + 48.6) / (-2.5)) * dfac * freq * wave_cm #we want to convert from Luminosity to int(lambda*Lum_lambda*dlambda)
    Q[ind] = lum / hc #rate of ionising photons in s^-1.
  
    return Q


def plot_halpha_lf_evo(plt, outdir, obsdir, h0, LFs_halpha):

    volcorr = 3.0*np.log10(h0)
    xlf_obs  = xlf
 
    xtit="$\\rm log_{10}(L_{H\\alpha}/erg\\, s^{-1})$"
    ytit="$\\rm log_{10}(\Phi/{\\rm dex^{-1}} {\\rm Mpc}^{-3})$"

    xmin, xmax, ymin, ymax = 39, 45, -6, -1
    xleg = xmax - 0.3 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)

    fig = plt.figure(figsize=(5,14))

    
    subplots = (511, 512, 513, 514, 515)
    idx = (0, 1, 2, 3, 4)
    zs  = (0, 1, 2, 3, 4)
    labels= ('z=0', 'z=0.25', 'z=1', 'z=1.5', 'z=2')
  
    corrm_obs = np.log10(h0/0.7)
    corry_obs = 3.0*np.log10(h0/0.7)

    for subplot, idx, z in zip(subplots, idx, zs):

        ax = fig.add_subplot(subplot)
        ytitplot = ytit
        if (idx == 4):
            xtitplot = xtit
        else:
            xtitplot = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitplot, ytitplot, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labels[idx])

        if(idx == 0):
           file = obsdir+'/lf/Halpha/GAMA_Halpha.dat'
           lm,p,dpl,dph = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           yobs = p[0:15]
           ydn  = p[0:15]-dpl[0:15]
           yup  = p[0:15]+dph[0:15]
           ax.errorbar(lm[0:15]+7+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Gunawardhana+2013")

        if(idx == 1):
           file = obsdir+'/lf/Halpha/GAMA_Halpha.dat'
           lm,p,dpl,dph = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           print(p.shape)
           yobs = p[24:43]
           ydn  = p[24:43]-dpl[24:43]
           yup  = p[24:43]+dph[24:43]
           ax.errorbar(lm[24:43]+7+corrm_obs, yobs+corry_obs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

           file = obsdir+'/lf/Halpha/Ly_0.24.dat'
           lm,p,dp = np.loadtxt(file,usecols=[0, 1, 2],unpack=True)
           yobs = np.log10(p)
           ydn  = np.log10(p-dp)
           yup  = np.log10(p+dp)
           ax.errorbar(lm+corrm_obs, yobs+corry_obs, yerr=[abs(yobs-ydn),abs(yup-yobs)], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='^',label="Ly+2007")

        if(idx == 2):
           file = obsdir+'/lf/Halpha/Sobral13_z0.8.dat'
           lm,p,dp = np.loadtxt(file,usecols=[0,5,6],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=dp, ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s', label='Sobral+2013')
           file = obsdir+'/lf/Halpha/Ly_2010_z0p8.dat'
           lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
           ax.errorbar(lm+corrm_obs, p-0.3+corry_obs, yerr=abs(p-dp), ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*', label='Ly+2010')

        if(idx == 3):
           file = obsdir+'/lf/Halpha/Sobral13_z1.4.dat'
           lm,p,dp = np.loadtxt(file,usecols=[0,5,6],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=dp, ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')

        if(idx == 4):
           file = obsdir+'/lf/Halpha/Sobral13_z2.23.dat'
           lm,p,dp = np.loadtxt(file,usecols=[0,5,6],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=dp, ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')
           file = obsdir+'/lf/Halpha/Haynes_2010_z2p2.dat'
           lm,p,dpl,dph = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
           ax.errorbar(lm+corrm_obs, p+corry_obs, yerr=[p-dpl, dph-p], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='p', label = 'Hayes+2010')

        #Predicted LF
        ind = np.where(LFs_halpha[z,0,:] < 0.)
        y = LFs_halpha[z,0,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'k', linewidth=3, label = 'total' if idx == 0 else None)
        ind = np.where(LFs_halpha[z,1,:] < 0.)
        y = LFs_halpha[z,1,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'b', linewidth=3, label = 'disks' if idx == 0 else None)
        ind = np.where(LFs_halpha[z,2,:] < 0.)
        y = LFs_halpha[z,2,ind]+volcorr-np.log10(dm)
        ax.plot(xlf_obs[ind],y[0],'r', linewidth=3, label = 'bulges' if idx == 0 else None)

        if(idx == 0):
            common.prepare_legend(ax, ['k', 'b', 'r', 'grey','grey','grey'], loc=3)
        else:
            common.prepare_legend(ax, ['grey','grey','grey'], loc=3)
        #ind = np.where(LFs_nodust[z,4,band,:] < 0.)
        #y = LFs_nodust[z,4,band,ind]+volcorr-np.log10(dm)
        #ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

    plt.tight_layout()
    common.savefig(outdir, fig, "Halpha_luminosity_function_evolution.pdf")

def prepare_data(hdf5_data, phot_data_nod, LFs_Halpha, index, nbands):
   
    #star_formation_histories and SharkSED have the same number of galaxies in the same order, and so we can safely assume that to be the case.
    #to select the same galaxies in galaxies.hdf5 we need to ask for all of those that have a stellar mass > 0, and then assume that they are in the same order.

    (h0, volh, mdisk, mbulge, mhalo, mshalo, typeg, age,
     sfr_disk, sfr_burst, id_gal) = hdf5_data

    #components:
    #(len(my_data), 2, 2, 5, nbands)
    #0: disk instability bulge
    #1: galaxy merger bulge
    #2: total bulge
    #3: disk
    #4: total
    seds_bulge_dummy = phot_data_nod[0]
    ngals = len(seds_bulge_dummy[0,:])
    SEDs_nodust = np.zeros(shape = (5,nbands,ngals))
    SEDs_nodust[0,:] = phot_data_nod[0]
    SEDs_nodust[1,:] = phot_data_nod[1]
    SEDs_nodust[2,:] = phot_data_nod[2]
    SEDs_nodust[3,:] = phot_data_nod[3]
    SEDs_nodust[4,:] = phot_data_nod[4]

    ion_mag = SEDs_nodust[4,27,:] #27 == Band_ionising_photons
    Lhalpha = 1.37e-12 * ionising_photons(ion_mag, 912.0) #in erg/s
    ion_mag_d = SEDs_nodust[3,27,:] #27 == Band_ionising_photons
    Lhalpha_d = 1.37e-12 * ionising_photons(ion_mag_d, 912.0) #in erg/s
    ion_mag_b = SEDs_nodust[2,27,:] #27 == Band_ionising_photons
    Lhalpha_b = 1.37e-12 * ionising_photons(ion_mag_b, 912.0) #in erg/s

    print(Lhalpha)

    ind = np.where(Lhalpha > 1e35)
    H, bins_edges = np.histogram(np.log10(Lhalpha[ind]),bins=np.append(mbins,mupp))
    LFs_Halpha[index,0,:] = H

    ind = np.where(Lhalpha_d > 1e35)
    H, bins_edges = np.histogram(np.log10(Lhalpha_d[ind]),bins=np.append(mbins,mupp))
    LFs_Halpha[index,1,:] = H

    ind = np.where(Lhalpha_b > 1e35)
    H, bins_edges = np.histogram(np.log10(Lhalpha_b[ind]),bins=np.append(mbins,mupp))
    LFs_Halpha[index,2,:] = H


def main(model_dir, outdir, redshift_table, subvols, obsdir):

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()

    Variable_Ext = True

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age',
                           'sfr_disk', 'sfr_burst', 'id_galaxy')}

    fields_sed_nod = {'SED/ab_nodust': ('bulge_d','bulge_m','bulge_t','disk','total')}

    z = (0, 0.25, 1, 1.5, 2.0)  
    snapshots = redshift_table[z]

    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5" 
    # Create histogram
    for index, snapshot in enumerate(snapshots):
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)

        if(Variable_Ext == False):
           seds_nod = common.read_photometry_data(model_dir, snapshot, fields_sed_nod, subvols)
        else:
           seds_nod = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_nod, subvols, file_hdf5_sed)

        nbands = len(seds_nod[0]) 

        if(index == 0):
            LFs_halpha   = np.zeros(shape = (len(z), 3, len(mbins)))

        prepare_data(hdf5_data, seds_nod, LFs_halpha, index, nbands)

        h0, volh = hdf5_data[0], hdf5_data[1]
        if(volh > 0.):
            LFs_halpha[index,:]   = LFs_halpha[index,:]/volh

    # Take logs
    ind = np.where(LFs_halpha > 0.)
    LFs_halpha[ind] = np.log10(LFs_halpha[ind])

    if(Variable_Ext):
       outdir = os.path.join(outdir, 'eagle-rr14')

    plot_halpha_lf_evo(plt, outdir, obsdir, h0, LFs_halpha)

if __name__ == '__main__':
    main(*common.parse_args())
