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
"""Size plots"""

import functools

import numpy as np

import common
import utilities_statistics as us


zlist=np.array([3.7248038025221, 3.95972, 4.465197621546, 4.73693842543988, 5.02220991014863, 5.2202206934302, 5.52950356184419, 5.74417977285603, 5.96593, 6.19496927748119, 6.55269895697227, 7.05756323172746, 7.45816170313544, 7.73629493731708, 8.02352,8.32018565809831, 8.47220854014322, 8.78358705435761, 8.94312532315157, 9.27010372804765, 9.437541750167, 9.78074128377067, 9.95655])

##################################
#Constants
RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun
c_light  = 299792458.0 #m/s
PI       = 3.141592654

mlow = 6.5
mupp = 12.5
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

dmobs = 0.4
mbins_obs = np.arange(mlow,mupp,dmobs)
xmf_obs = mbins_obs + dmobs/2.0

vlow = 1.0
vupp = 3.0
dv   = 0.1
vbins = np.arange(vlow,vupp,dv)
xv    = vbins + dv/2.0

btbins = [0, 0.05, 0.95, 1]

btlow = 0.0
btupp = 1.0
dbt   = 0.05
btbins2 = np.arange(btlow,btupp,dbt)
xbt    = btbins2 + dbt/2.0

def smooth(x, y, ndeg):
    fit = np.polyfit(x,y,ndeg)
    print(fit) 
    y_smooth = np.zeros(shape = len(x))
    for j in range(0,ndeg+1):
        y_smooth[:] = y_smooth[:] + fit[j] * x[:]**(ndeg -j)

    return y_smooth

def prepare_data(hdf5_data, seds, index, zsnap, obsdir):

    (h0, volh, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     sfr_disk, sfr_bulge, typeg,  mgas_disk, mgas_bulge, matom_disk, mmol_disk, matom_bulge, mmol_bulge, mvir_hosthalo, idtree) = hdf5_data

    SEDs_dust_total = seds[4] 
    sfr_gals = (sfr_disk + sfr_bulge) / 1e9 / h0

    def calculate_lir(SEDs_dust_total, sfrs, obsdir):
        file = obsdir+'/Models/Shark_SED_bands.dat'
        lambda_bands = np.loadtxt(file,usecols=[0],unpack=True)
        freq_bands   = c_light / (lambda_bands * 1e-10) #in Hz

        print(SEDs_dust_total.shape)
        bands = [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 30] #from 8mu to 1000mu
        lir = np.zeros(shape = (len(sfrs)))
        fluxes = 10.0**(SEDs_dust_total[:,0,:] / -2.5) * 3631.0 * 1e-23 * 4.0 * PI * (3.086e19)**2.0 / 3.846e33 #in Lsun Hz−1
        print(fluxes.shape) 
        for i in range(0,len(sfrs)):
            for b in range(0,len(bands)-1):
                delta_freq = abs(freq_bands[bands[b]] - freq_bands[bands[b+1]])
                lir[i] = lir[i] + fluxes[bands[b], i] * delta_freq + abs(fluxes[bands[b+1], i] - fluxes[bands[b], i]) * delta_freq
        return lir
    ind = np.where(sfr_gals > 50)
    if(len(sfr_gals[ind]) > 0):
       lir_50sfr = calculate_lir(SEDs_dust_total[:,ind], sfr_gals[ind], obsdir)
       ind = np.where(lir_50sfr > 3e12)
       numberdensity = (len(lir_50sfr[ind]) + 0.0) / (volh * h0**3.0)
       print('number density of galaxies with LIR>3e12Msun is', numberdensity, ' at redshift ', zsnap)

    return(volh, h0)
    
def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'sfr_disk', 'sfr_burst', 'type', 
                           'mgas_disk', 'mgas_bulge','matom_disk', 'mmol_disk', 
                           'matom_bulge', 'mmol_bulge', 'mvir_hosthalo', 'id_halo_tree')}

    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5"
    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}
    fields_sed_nod = {'SED/ab_nodust': ('disk','total')}
    fields_sed_ap = {'SED/ap_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}

    #bands of interest band-7, band-6, band-4
    selec_alma = (29, 30, 32)
    flux_selec = (1e-10, 1e-2, 1e-1, 1.0) #to look at 0.01<S<0.1, 0.1<S<1, S>1mJy

    for index, snapshot in enumerate(redshift_table[zlist]):
        print("Will read snapshot %s" % (str(snapshot)))
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        seds = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed, subvols, file_hdf5_sed)

        (volh, h0) = prepare_data(hdf5_data, seds, index, zlist[index], obsdir)

if __name__ == '__main__':
    main(*common.parse_args())
