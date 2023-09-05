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
from scipy import interpolate
import common
import os
import math 
import h5py

##################################

# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3
c_light = 299792458.0 #m/s
Sigma0=1e-10 #yr^{-1}*zsun
FX0=1e-2 #erg/s/cm^2
Lsun = 3.839e-7 #in 1e40 erg/s
sigma_gas = 20.0 #km/s for CO

thresh_thin_disk = 0.01
thresh_super_edd = 1.0

gammaSFR = 1.0
alphaFx  = 1.0 
Av = 4

zsun = 0.0189

def prepare_data(hdf5_data, index, model_dir, snapshot, subvol, obsdir):

    #read ascii table from Bayet et al. (2011)
    datanu = np.genfromtxt(os.path.join(obsdir, 'Models', 'CO','emlines_carbonmonoxide_restnu.data'), delimiter=' ', skip_header=11)
    nuco = datanu

    dataPDR  = np.genfromtxt(os.path.join(obsdir, 'Models', 'CO','emlines_carbonmonoxide_Bayet11.data'), delimiter=' ', skip_header=10)

    CRs = dataPDR[:,5]/5E-17 #normalize cosmic rays rate by 5e-17 s^{-1}.
    Xconv_Av3 = dataPDR[:,6:16]/1e20 #normalize conversion factors by 10^20 cm^(-2) (K km/s)^(-1)
    Xconv_Av8 = dataPDR[:,16:26]/1e20 #normalize conversion factors by 10^20 cm^(-2) (K km/s)^(-1)

    #Convert all the relevant quantities to logarithm as the 3D interpolation is done in the logarithmic space.
    Xconv_Av3 = np.log10(Xconv_Av3)
    Xconv_Av8 = np.log10(Xconv_Av8)
 
    Zmod = np.log10(dataPDR[:,4])
    GUV  = np.log10(dataPDR[:,0])
    CRs  = np.log10(CRs)

    nh   = dataPDR[:,3]
    ind = np.where(nh == 1e4)
    interpolator = interpolate.LinearNDInterpolator(list(zip(Zmod[ind], GUV[ind], CRs[ind])), np.squeeze(Xconv_Av3[ind]))
    interpolator_nn = interpolate.NearestNDInterpolator(list(zip(Zmod[ind], GUV[ind], CRs[ind])), np.squeeze(Xconv_Av3[ind]))

    MinZ = min(Zmod) #define minimum metallicity probed by the models.
    MaxZ = max(Zmod) #define maximum metallicity probed by the models.

    MinGUV = np.min(GUV)  #define minimum GUV probed by the models.
    MaxGUV = np.max(GUV)  #define maximum GUV probed by the models.

    MinCRs = np.min(CRs) #define minimum CRs probed by the models.
    MaxCRs = np.max(CRs) #define maximum CRs probed by the models.

    # read galaxy information in lightcone
    (h0, _, idgal, mmol_b, mmol_d, rd, rb, mzd, mzb, sfr_d, sfr_b, mgd, mgb, 
    mbh, mbh_acc_hh, mbh_acc_sb, mdisk, mbulge) = hdf5_data

    #define metallicities for disk and bulge
    zd = np.zeros(shape = len(mzd))
    zb = np.zeros(shape = len(mzb))
    ind = np.where(mgd > 0)
    zd[ind] = mzd[ind]/mgd[ind]
    ind = np.where(mgb > 0)
    zb[ind] = mzb[ind]/mgb[ind]

    # The calculation below is done to ultimately compute the amout of X-ray flux in the galaxy:

    # Eddington luminosity calculation
    Ledd = 1.26e6 * (mbh/h0/1e8) #in 1e40 ergs/s

    # Eddington accretion rate from Eq 8 in Griffin et al. (208)
    macc_edd = Ledd/(0.1*pow(c_light*1e2,2.0)) * 1.577e+23 #in units of Msun/Gyr

    mnorm = (mbh_acc_hh + mbh_acc_sb)/h0/macc_edd #accretion rate normalized by Eddington rate

    Lbol = np.zeros(len(Ledd))
    Lthin_disk = np.zeros(len(Ledd))

    # define bolometric luminosities using Eqs 1 in Amarantidis et al. (2019)
    ind = np.where(mbh_acc_hh+mbh_acc_sb > 0) 
    Lthin_disk[ind] = 0.1*pow(c_light*1e2, 2.0) *(mbh_acc_hh[ind]+mbh_acc_sb[ind])/h0 * 6.329113924050633e-24 #in 1e40 ergs/s

    # thin disks
    ind = np.where((mnorm > thresh_thin_disk) & (mnorm < thresh_super_edd))
    Lbol[ind] = Lthin_disk[ind]
    # super-eddington
    ind = np.where(mnorm > thresh_super_edd)
    Lbol[ind] = thresh_super_edd *  (1.0 + np.log(mnorm[ind]/thresh_super_edd)) * Ledd[ind] #in 1e40 ergs/s
    # ADAfs adopting alpa_ADAF=0.1 and beta = 0.81 and r_lso = 6 (spin =1)
    ind = np.where((mnorm < thresh_thin_disk) & (mnorm > 0))
    Lbol[ind] = 0.2 * Lthin_disk[ind] * (mnorm[ind]/1e-2) * 1.63 
 
    # L-hard-xrays Eq 34 from Griffin et al. (2018)
    Lrat = np.log10(Lbol/Lsun)  
    Lx = pow(10.0, -1.54 - 0.24*Lrat - 0.012 * pow(Lrat, 2.0) + 0.0015 * pow(Lrat, 3.0)) * Lbol #in 1e40 erg/s

    ind = np.where(Lbol <= 0)
    Lx[ind] = 0

    # define relevant quantities we need for the calculation of CO
    SFRtot = (sfr_d + sfr_b)/1e9/h0 #Msun/yr
    SFRburst = sfr_b/1e9/h0 #Msun/yr
    SFRdisk = sfr_d/1e9/h0 #Msun/yr

    r50_disk = rd * 1e3/ h0 #kpc
    r50_bulge = rb * 1e3/ h0 #kpc

    zcoldg_d = np.log10(zd/zsun)
    z_zero = np.where(zd <= 0)
    zcoldg_d[z_zero] = MinZ
    zcoldg_d = np.clip(zcoldg_d, MinZ, MaxZ)

    zcoldg_b = np.log10(zb/zsun)
    ind = np.where(zb <= 0)
    zcoldg_b[ind] = MinZ
    zcoldg_b = np.clip(zcoldg_b, MinZ, MaxZ)
 
    # calculation of quantities that go directly into the CO computation
    # calculation UV radiation field.
    def get_guv(sfr, mgas, z):
        guv = (sfr / mgas / (z/zsun)) / Sigma0
        guv = gammaSFR * np.log10(guv)    
        is_zero = np.where((mgas <= 0) | (sfr <= 0))
        guv[is_zero] = MinGUV
        return np.clip(guv, MinGUV, MaxGUV)

    guv_disk = get_guv(SFRdisk, mgd/h0, zd)
    guv_bulge = get_guv(SFRburst, mgb/h0, zb)

    # calculation X-ray radiation field.
    def get_xray(Lx,r):
        
        xray_field = np.zeros(len(Lx))
        ind = np.where((Lx > 0) & (Lx < 1e10) & (r > 0))
        xray_field[ind] = Lx[ind] / (4.0 * PI * pow(r[ind],2.0)) * 0.0010500455929796473 #in erg/s/cm^2
        xray_field[ind] = np.log10(xray_field[ind]/ FX0) #in solar units
        return np.clip(xray_field, MinCRs, MaxCRs) 

    xray_disk = np.zeros(len(Lx)) # assume no X-ray boost in disks
    xray_bulge = get_xray(Lx, rb/h0)

    # calculate CO emission of galaxies given some inputs
    def get_co_emissions(mmol, zcoldg, guv, fx):

        shape = len(mmol), 10
        CRRayFlux = fx #np.full((len(guv)), 1.)

        ind = np.where(mmol > 0)
        mcold = np.zeros(shape)
        for i in range(0,10):
            mcold[ind, i] = mmol[ind]

	# Interpolate linearly first
        # If extrapolation is needed we use the nearest neighbour interpolator
        xco = np.zeros(shape)
        xco[ind, :] = 10.0 ** interpolator(list(zip(zcoldg[ind], guv[ind], CRRayFlux[ind])))
        isnan = np.where(np.isnan(xco[:, 0]))
        xco[isnan, :] = 10.0 ** interpolator_nn(list(zip(zcoldg[isnan], guv[isnan], CRRayFlux[isnan])))

        lco = np.zeros(shape)
        lco[ind, :] = mcold[ind] * XH / 313./ xco[ind]
        for i in range(0, 10):
            lco[ind, i] = lco[ind, i] * pow(i + 1.0, 2.0)
        return lco

    # calculate CO luminosity coming from the disk and the bulge in [Jy km/s Mpc^2]
    LCOb = get_co_emissions(mmol_b/h0, zcoldg_b, guv_bulge, xray_bulge)
    LCOd = get_co_emissions(mmol_d/h0, zcoldg_d, guv_disk, xray_disk)
    # will write the hdf5 files with the CO SLEDs and relevant quantities
    # will only write galaxies with mstar>0 as those are the ones being written in SFH.hdf5
    ind = np.where( (mdisk +  mbulge) > 0)

    # will write the hdf5 files with the CO SLEDs and relevant quantities
    file_to_write = os.path.join(model_dir, str(snapshot), str(subvol), 'CO_SLED.hdf5')
    print ('Will write extinction to %s' % file_to_write)
    hf = h5py.File(file_to_write, 'w')
    hf.create_dataset('galaxies/id_galaxy', data=idgal[ind])
    hf.create_dataset('galaxies/LCO_disk', data=LCOd[ind])
    hf.create_dataset('galaxies/LCO_bulge', data=LCOb[ind])
    hf.create_dataset('galaxies/Lum_AGN_HardXray', data=Lx[ind])
    hf.create_dataset('frequency_co_rest', data=nuco)
    hf.close()

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):

    snapshots = range(30,200)
    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()

    fields = {'galaxies': ('id_galaxy', 'mmol_bulge', 'mmol_disk', 'rgas_disk',
                           'rgas_bulge', 'mgas_metals_disk', 'mgas_metals_bulge', 'sfr_disk', 'sfr_burst',
                           'mgas_disk', 'mgas_bulge','m_bh','bh_accretion_rate_hh','bh_accretion_rate_sb',
                           'mstars_disk', 'mstars_bulge')}

    for index,snapshot in enumerate(snapshots): 
       for subv in subvols:
           hdf5_data = common.read_data(model_dir, snapshot, fields, [subv])
           prepare_data(hdf5_data, index, model_dir, snapshot, subv, obs_dir)

if __name__ == '__main__':
    main(*common.parse_args())
