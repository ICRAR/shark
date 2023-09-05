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
mlow = 6.0
mupp = 12.0
dm = 0.2
mbins = np.arange(mlow, mupp, dm)
xmf = mbins + dm/2.0

zsun = 0.0189
MpctoKpc = 1e3

#model of Mattson et al. (2014) for the dependence of the dust-to-metal mass ratio and metallicity X/H.
corrfactor_dm = 2.0
polyfit_dm = [ 0.00544948, 0.00356938, -0.07893235,  0.05204814,  0.49353238] 

#choose dust model between mm14, rr14 and constdust
m14 = False
rr14 = True
constdust = False
rr14xcoc = False
random_perturbation = True

#read EAGLE tables
sdust_eaglet, taumed_eagle, taulow_eagle, tauhigh_eagle = common.load_observation('../data', 'Models/EAGLE/Tau5500-Trayford-EAGLE.dat', [0,1,2,3])
sdust_eaglem, mmed_eagle, mlow_eagle, mhigh_eagle = common.load_observation('../data/','Models/EAGLE/CFPowerLaw-Trayford-EAGLE.dat', [0,1,2,3])

def interp (sdust_eagle, med_eagle, low_eagle, high_eagle):

    med = np.zeros(shape = (2,len(sdust_eagle)-1))
    low = np.zeros(shape = (2,len(sdust_eagle)-1))
    hig = np.zeros(shape = (2,len(sdust_eagle)-1))
   
    for i in range(0,len(sdust_eagle)-1):
        delta_dust = sdust_eagle[i+1] - sdust_eagle[i]
        med[0,i] = (med_eagle[i+1] - med_eagle[i] ) / delta_dust
        med[1,i] =  med_eagle[i+1]- med[0,i] * sdust_eagle[i+1]
        low[0,i] = (low_eagle[i+1] - low_eagle[i]) / delta_dust
        low[1,i] =  low_eagle[i+1]- low[0,i] * sdust_eagle[i+1]
        hig[0,i] = (high_eagle[i+1] - high_eagle[i]) / delta_dust
        hig[1,i] =  high_eagle[i+1]- hig[0,i] * sdust_eagle[i+1]

    return (med, low, hig)

(m_med, m_low, m_hig) = interp (sdust_eaglet, taumed_eagle, taulow_eagle, tauhigh_eagle)
(s_med, s_low, s_hig) = interp (sdust_eaglem, mmed_eagle, mlow_eagle, mhigh_eagle)

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

# define tau diffuse
def tau_diff (md, rd, hd, h0):

    tau    = np.zeros(shape = len(md))
    sigma  = np.zeros(shape = len(md))

    ind = np.where((md > 0) & (rd > 0) & (hd > 0))
    sigma[ind] = np.log10(md[ind]/h0 / (2.0 * 3.1416 * rd[ind]*MpctoKpc/h0 * hd[ind]*MpctoKpc/h0)) #in Msun/kpc^2
    # cap surface density of dust to physical values based on EAGLE
    sigma[ind] = np.clip(sigma[ind],0, 12)

    #interpolate linearly in dust surface density going through the list of values in the EAGLE table
    for i in range(0,len(sdust_eaglet)-1):
        if(i == 0):
           selecinrage = np.where(sigma <= sdust_eaglet[i+1])
        elif(i == len(sdust_eaglet)-2):
           selecinrage = np.where(sigma >= sdust_eaglet[i])
        else:
           selecinrage = np.where((sigma >= sdust_eaglet[i]) & (sigma < sdust_eaglet[i+1]))

        tau[selecinrage] = m_med[0,i] * sigma[selecinrage] + m_med[1,i]
        tlow = abs(tau[selecinrage] - (m_low[0,i] * sigma[selecinrage] + m_low[1,i]))
        thigh= abs( (m_hig[0,i] * sigma[selecinrage] + m_hig[1,i]) - tau[selecinrage])
        var_gauss = (tlow + thigh) * 0.5
        pert = np.random.randn(len(var_gauss)) * np.sqrt(var_gauss)
        if random_perturbation == True:
           tau[selecinrage] = tau[selecinrage] + pert

    # cap it to maximum and minimum values in EAGLE
    tau = np.clip(tau, 1e-6, 5) 

    return (tau, sigma) 

# define clump tau

def tau_clump(mz,mg,h0, sigmag, tau_diff):
    sigmaclump = np.zeros(shape = len(mg))
    sigmaclump[:] = 85.0*1e6 #in Msun/kpc^3
    ind = np.where(sigmaclump < sigmag)
    sigmaclump[ind] = sigmag[ind]
    tau = np.zeros(shape = len(mz))
    ind = np.where((mz > 0) & (mg > 0))
    (md, DToM_MW)  = dust_mass(mz[ind],mg[ind],h0)
    norm = 85.0 * 1e6 * DToM_MW * zsun #dust surface density of clumps
    tau[ind] = 1.0 * (sigmaclump[ind] * md/(mg[ind]/h0) / norm)
    ind = np.where(tau < tau_diff)
    tau[ind] = tau_diff[ind]
    # cap it to maximum and minimum values in EAGLE but also forcing the clump tau to be at least as high as the diffuse tau
    tau = np.clip(tau, 1e-6, 5)
    return tau


def tau_clump2(mz,mg,h0):
    tau = np.zeros(shape = len(mz))
    ind = np.where((mz > 0) & (mg > 0))
    (md, DToM_MW)  = dust_mass(mz[ind],mg[ind],h0)
    tau[ind] = 1.5 * (md/mz[ind]/DToM_MW)
    # cap it to maximum and minimum values in EAGLE but also forcing the clump tau to be at least as high as the diffuse tau
    tau = np.clip(tau, 1e-6, 5)
    return tau

# define slope diffuse
def slope_diff (md, rd, hd, h0):

    m    = np.zeros(shape = len(md))
    sigma  = np.zeros(shape = len(md))

    ind = np.where((md > 0) & (rd > 0) & (hd > 0))
    sigma[ind] = np.log10(md[ind]/h0 / (2.0 * 3.1416 * rd[ind]*MpctoKpc/h0 * hd[ind]*MpctoKpc/h0)) #in Msun/kpc^2
    # cap surface density of dust to physical values based on EAGLE
    sigma[ind] = np.clip(sigma[ind],0, 12)

    #interpolate linearly in dust surface density going through the list of values in the EAGLE table
    for i in range(0,len(sdust_eaglem)-1):
        if(i == 0):
           selecinrage = np.where(sigma <= sdust_eaglem[i+1])
        elif(i == len(sdust_eaglet)-2):
           selecinrage = np.where(sigma >= sdust_eaglem[i])
        else:
           selecinrage = np.where((sigma >= sdust_eaglem[i]) & (sigma < sdust_eaglem[i+1]))

        m[selecinrage] = s_med[0,i] * sigma[selecinrage] + s_med[1,i]
        tlow = abs(m[selecinrage] - (s_low[0,i] * sigma[selecinrage] + s_low[1,i]))
        thigh= abs((s_hig[0,i] * sigma[selecinrage] + s_hig[1,i]) - m[selecinrage])
        var_gauss = (tlow + thigh) * 0.5
        pert = np.random.randn(len(var_gauss)) * np.sqrt(var_gauss)
        if random_perturbation == True:
           m[selecinrage] = m[selecinrage] + pert

    # cap it to maximum and minimum values in EAGLE
    m = np.clip(m, -3, -0.001)
    return m 


def prepare_data(hdf5_data, index, model_dir, snapshot, subvol):

    bin_it = functools.partial(us.wmedians, xbins=xmf)

    # Unpack data
    (h0, _, typeg, rgasd, rgasb, mHId, mH2d, mgasd, mHIb, mH2b, mgasb, mzd, mzb, mdisk, mbulge, sfrd, sfrb, idgal) = hdf5_data
    XH = 0.72
    h0log = np.log10(float(h0))

    sigma_g_d = mgasd/h0/(2.0 * 3.1416 * (rgasd/h0*1e3)**2.0)
    sigma_g_b = mgasb/h0/(2.0 * 3.1416 * (rgasb/h0*1e3)**2.0)


    (mdustd, DToM_MW) = dust_mass(mzd, mgasd, h0)
    (mdustb, DToM_MW) = dust_mass(mzb, mgasb, h0)

    #random numbers from 0 to 1
    sin_inclination = np.random.rand(len(mdustd)) 
    inclination = np.arcsin(sin_inclination) * 180.0/3.1416 #degrees
    bd  = sin_inclination*(rgasd - rgasd/7.3)  +  rgasd/7.3 #scaleheight at r50

    (tau_dust_bulge, sigmab) = tau_diff(mdustb, rgasb, rgasb, h0)
    (tau_dust_disk, sigmad) = tau_diff(mdustd, rgasd, bd, h0)
    
    tau_clump_bulge = tau_clump(mzb, mgasb, h0, sigma_g_b, tau_dust_bulge)
    tau_clump_disk  = tau_clump(mzd, mgasd, h0, sigma_g_d, tau_dust_disk)
    slope_dust_bulge = slope_diff(mdustb, rgasb, rgasb, h0) 
    slope_dust_disk  = slope_diff(mdustd, rgasd, bd, h0)

    tau_clump_bulge2 = tau_clump2(mzb, mgasb, h0)
    tau_clump_disk2  = tau_clump2(mzd, mgasd, h0)


    # will write the hdf5 files with the CO SLEDs and relevant quantities
    # will only write galaxies with mstar>0 as those are the ones being written in SFH.hdf5
    ind = np.where( (mdisk +  mbulge) > 0)
    file_to_write = os.path.join(model_dir, str(snapshot), str(subvol), 'extinction-eagle-rr14.hdf5')
    print ('Will write extinction to %s' % file_to_write)
    hf = h5py.File(file_to_write, 'w')
    
    hf.create_dataset('galaxies/tau_screen_disk', data=tau_dust_disk[ind])
    hf.create_dataset('galaxies/tau_screen_bulge', data=tau_dust_bulge[ind])
    hf.create_dataset('galaxies/tau_birth_disk', data=tau_clump_disk[ind])
    hf.create_dataset('galaxies/tau_birth_bulge', data=tau_clump_bulge[ind])
    hf.create_dataset('galaxies/pow_screen_disk', data=slope_dust_disk[ind])
    hf.create_dataset('galaxies/pow_screen_bulge', data=slope_dust_bulge[ind])
    hf.create_dataset('galaxies/id_galaxy', data=idgal[ind])
    hf.create_dataset('galaxies/inclination', data=inclination[ind])
    hf.close()

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):

    #zlist = np.arange(2,10,0.25)
    #zlist = (0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 0, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 9, 10)

    zlist_given = False
    if(zlist_given):
        zlist = [0.381963715160695, 1.77053590476006]
    else:
        snap_list = [269, 224, 213, 205, 188, 174, 153, 140, 129, 120, 111, 104, 91, 82, 75, 69, 63, 58]
        #[199, 185, 179, 174, 164, 156, 149, 142, 136, 131, 113, 100, 88, 79, 70, 63, 57, 51]

    #269 224 213 205 188 174 153 140 129 120 111 104 91 82 75 69 63 58
    #zlist = [0, 0.254144, 0.5, 0.75744098, 1.00678003, 1.49550998, 2.00202990, 
    #zlist = [0.381963715160695, 1.77053590476006] 
            #0, 0.24944700, 0.38672200, 0.49594200, 0.75744098, 1.00678003, 1.49550998, 2.00202990, 2.51012993, 2.98918009, 3.53362989, 4.00793982, 5.02441978, 6.01074982, 7.00544024, 7.96958017, 9.04985046, 10.04880047] 
    #1.0, 1.5, 2.0, 3.0, 3.95, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0] #0.25, 0.38, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0]
            #0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 3.95, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0] #0.25, 0.38, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0]
    #8.02352, 8.94312532315157, 9.95655, 10.5013916683919, 12.520639255824]
    #0, 0.194739, 0.254144, 0.359789, 0.450678, 0.8, 0.849027, 0.9, 1.20911, 1.28174, 1.39519, 1.59696, 2.00392, 2.47464723643932, 2.76734390952347, 3.01916, 3.21899984389701, 3.50099697082904, 3.7248038025221, 3.95972, 4.465197621546, 4.73693842543988, 5.02220991014863, 5.52950356184419, 5.96593, 6.55269895697227, 7.05756323172746, 7.45816170313544, 8.02352, 8.94312532315157, 9.95655, 10.5013916683919, 12.520639255824] 
    #0.033] 
    #0, 0.25, 0.38, 0.49, 0.76, 1.00, 1.245, 1.51, 1.77, 2.00, 3.01, 3.95, 5.0, 6.0, 7.0, 8.0, 9.0, 10]
    #[0,0.194739, 0.254144, 0.359789, 0.450678, 0.8, 0.849027, 0.9, 1.20911, 1.28174, 1.39519, 1.59696, 2.00392, 2.47464723643932, 2.76734390952347, 3.01916, 3.21899984389701, 3.50099697082904, 3.7248038025221, 3.95972, 4.465197621546, 4.73693842543988, 5.02220991014863, 5.52950356184419, 5.96593, 6.55269895697227, 7.05756323172746, 7.45816170313544, 8.02352, 8.94312532315157, 9.95655]
    #[0.016306640039433, 0.066839636933135, 0.084236502339783, 0.119886040396529, 0.138147164704691, 0.175568857770275, 0.214221447279112, 0.23402097095238, 0.274594901875312, 0.316503156974571]

    plt = common.load_matplotlib()
    fields = {'galaxies': ('type', 'rgas_disk', 'rgas_bulge', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge', 'mgas_metals_disk', 
                           'mgas_metals_bulge', 'mstars_disk', 'mstars_bulge','sfr_disk','sfr_burst','id_galaxy')}

    if(zlist_given): 
       for index, snapshot in enumerate(redshift_table[zlist]):
           for subv in subvols:
               hdf5_data = common.read_data(model_dir, snapshot, fields, [subv])
               prepare_data(hdf5_data, index, model_dir, snapshot, subv)
    else:
        for index, snapshot in enumerate(snap_list):
            for subv in subvols:
               hdf5_data = common.read_data(model_dir, snapshot, fields, [subv])
               prepare_data(hdf5_data, index, model_dir, snapshot, subv)

if __name__ == '__main__':
    main(*common.parse_args())
