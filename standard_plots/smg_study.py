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


zlist=np.array([0.194739, 0.254144, 0.359789, 0.450678, 0.8, 0.849027, 0.9, 1.0, 1.10447, 1.20911, 1.28174, 1.39519, 1.51429, 1.59696, 1.68234, 1.81572, 1.90829, 2.00392, 3.01916, 3.95972, 5.96593, 8.02352, 9.95655])

##################################
#Constants
RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun

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


def prepare_data(hdf5_data, seds, seds_nod, seds_ap, index, sfr_z, mmol_z, sfr_tot, mmol_tot, flux_selec, selec_alma):

    (h0, volh, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     sfr_disk, sfr_bulge, typeg,  mgas_disk, mgas_bulge, matom_disk, mmol_disk, matom_bulge, mmol_bulge) = hdf5_data
    
    sfr_tot[index] = np.sum((sfr_disk + sfr_bulge) / 1e9 / h0)
    mmol_tot[index] = np.sum((mmol_bulge + mmol_disk) / h0)

    mstars_tot = (mdisk+mbulge)/h0
    bt = mbulge / (mdisk+mbulge)
    ind = np.where(mstars_tot > 0)
    mstars = mstars_tot[ind]
    SEDs_dust_total = seds[4]
    SEDs_vodust_total = seds_nod[1]
    SEDs_app = seds_ap[4]

    print(SEDs_app.shape)
    ind = np.where(mstars_tot > 0)
    sfrs_gals = (sfr_disk[ind] + sfr_bulge[ind]) / 1e9 / h0

    #calculate total SFR of galaxies selected in different ALMA bands and different fluxes
    #np.zeros(shape = (len(selec_alma), len(flux_selec), len(zlist)))
    for b in range(0,len(selec_alma)):
        flux_gals_band = 10.0**(SEDs_app[selec_alma[b],:] / -2.5) * 3631.0 * 1e3 #in mJy
        print("maximum flux ", max(flux_gals_band))
        for f in range(0,len(flux_selec)):
            if (f < 3):
                ind = np.where((flux_gals_band > flux_selec[f]) & (flux_gals_band <= flux_selec[f+1]))
            else:
                ind = np.where((flux_gals_band > flux_selec[f]) & (flux_gals_band < 1e10))
            sfr_z[b,f,index] = np.sum(sfrs_gals[ind])

    return(volh, h0)
    
def plot_sfr_contribution(plt, outdir, obsdir, sfr_z, sfr_tot):

    fig = plt.figure(figsize=(12,4.5))
    ytit = "$\\rm log_{10} (\\rm \\rho_{\\rm SFR}/ M_{\\odot} yr^{-1} cMpc^{-3})$"
    xtit = "redshift"
    xmin, xmax, ymin, ymax = 0, 10, -5, 0
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = [131, 132, 133]
    bands = ['band-7', 'band-6', 'band-4']
    cols = ['DarkBlue','MediumTurquoise','YellowGreen', 'Crimson']
    labels = ['$\\rm <0.01\\rm mJy$', '$\\rm  0.01-0.1\\rm mJy$', '$\\rm  0.1-1\\rm mJy$', '$\\rm  >1\\rm mJy$']

    for b in range(0,len(bands)):
        ax = fig.add_subplot(subplots[b])
        ytitle = ytit
        if(b > 0):
           ytitle=" "
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle, locators=(0.1, 1, 0.1, 1))
        ax.text(xleg,yleg,bands[b],fontsize=12) 
        for i in range(0,len(labels)):
            inp = np.where(sfr_z[b,i,:] != 0)
            y = sfr_z[b,i,inp]
            ax.plot(zlist[inp], y[0], linestyle='solid',color=cols[i], label=labels[i])
        ax.plot(zlist, sfr_tot,  linestyle='solid',color='k')
        common.prepare_legend(ax, cols, loc=4)
    common.savefig(outdir, fig, 'SFR_evolution_SMG_contribution.pdf')

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'sfr_disk', 'sfr_burst', 'type', 
                           'mgas_disk', 'mgas_bulge','matom_disk', 'mmol_disk', 
                           'matom_bulge', 'mmol_bulge')}

    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5"
    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}
    fields_sed_nod = {'SED/ab_nodust': ('disk','total')}
    fields_sed_ap = {'SED/ap_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}

    #bands of interest band-7, band-6, band-4
    selec_alma = (29, 30, 32)
    flux_selec = (1e-10, 1e-2, 1e-1, 0.5) #to look at 0.01<S<0.1, 0.1<S<1, S>1mJy

    sfr_z = np.zeros(shape = (len(selec_alma), len(flux_selec), len(zlist)))
    mmol_z = np.zeros(shape = (len(selec_alma), len(flux_selec), len(zlist)))
    sfr_tot = np.zeros(shape = (len(zlist)))
    mmol_tot = np.zeros(shape = (len(zlist)))

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        seds = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed, subvols, file_hdf5_sed)
        seds_nod = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed_nod, subvols, file_hdf5_sed)
        seds_ap = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed_ap, subvols, file_hdf5_sed)

        (volh, h0) = prepare_data(hdf5_data, seds, seds_nod, seds_ap, index, sfr_z, mmol_z, sfr_tot, mmol_tot, flux_selec, selec_alma)

    def take_log(x,v,h):
        x = x / (v * h**3.0)
        ind = np.where(x > 0)
        x[ind] = np.log10(x[ind])
        return x

    sfr_z = take_log(sfr_z, volh, h0)
    sfr_tot = take_log(sfr_tot, volh, h0)
    mmol_tot = take_log(mmol_tot, volh, h0)

    print(sfr_z[2,0,:], sfr_z[2,1,:], sfr_z[2,2,:])
    plot_sfr_contribution(plt, outdir, obsdir, sfr_z, sfr_tot)

if __name__ == '__main__':
    main(*common.parse_args())
