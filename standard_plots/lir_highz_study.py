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


zlist=np.array([0, 0.909822023685613, 2, 3, 4, 5])

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
Lsunwatts = 3.846e26

def prepare_data(hdf5_data, seds, seds_bands, fields_sed_bc, index, LFs_dust, obsdir):

    (h0, volh, mdisk, mbulge, sfr_disk, sfr_burst) = hdf5_data

    ind = np.where(mdisk + mbulge > 0)
    mstar = mdisk[ind] + mbulge[ind]

    bin_it = functools.partial(us.wmedians, xbins=xmf)

    seds_total = seds_bands[2]

    lir_disk = seds[0]
    lir_bulge = seds[1]
    lir_total = seds[2] #total absolute magnitudes with dust
    lir_bc_cont = fields_sed_bc[2]
    d2 = 4.0 * PI * (10 * 3.086e18)**2.0 #cm^2
    lir_1p4GHz = 10**((seds_total[12,:] + 48.6) / (-2.5)) / 1e7 * d2 #in W/Hz
    lir_3GHz = 10**((seds_total[11,:] + 48.6) / (-2.5)) / 1e7 * d2 #in W/Hz


    qIR = np.log10(lir_total[0,:]*Lsunwatts/3.75e12) - np.log10(lir_1p4GHz)
    qIR_3GHz = np.log10(lir_total[0,:]*Lsunwatts/3.75e12) - np.log10(lir_3GHz)
    ind = np.where((qIR > 0) & (qIR < 10))
    print(max(qIR[ind]), min(qIR[ind]), min(qIR_3GHz[ind]), max(qIR_3GHz[ind]), np.median(qIR[ind]))
    mstarbins = [8.5,9.5,10.5,11.5]
    dbin = 0.5
    for ms in mstarbins:
        ind = np.where((np.log10(mstar) >= ms - dbin) & (np.log10(mstar) < ms + dbin) & (qIR > 0) & (qIR < 10))
        print("Median qIR for stellar masses,", ms, " is", np.median(qIR[ind]))

    for ms in mstarbins:
        ind = np.where((np.log10(lir_total[0,:]) >= ms - dbin) & (np.log10(lir_total[0,:]) < ms + dbin) & (qIR > 0) & (qIR < 10))
        print("Median qIR for LIR,", ms, " is", np.median(qIR[ind]))

    temp_bc = 57.0
    temp_diff = 22.0
    temp_eff =  temp_bc * lir_bc_cont + temp_diff * (1.0 - lir_bc_cont)

    temp_eff = temp_eff[0]
    ind = np.where(temp_eff > temp_bc-1)
    print("number of galaxies with T=T_bc", len(temp_eff[ind]))
    print("all galaxies", len(temp_eff))
    ind = np.where(temp_eff < temp_diff+1)
    print("number of galaxies with T=T_diff", len(temp_eff[ind]))

    print(lir_total.shape)
    lir_total = lir_total[0]
    ind = np.where((lir_total > 1e9) & (temp_eff > temp_bc-1))
    print("number of galaxies with T=T_bc", len(temp_eff[ind]))
    ind = np.where((lir_total > 1e9))
    print("number of all galaxies", len(temp_eff[ind]))
    ind = np.where((lir_total > 1e9) & (temp_eff < temp_diff+1))
    print("number of galaxies with T=T_diff", len(temp_eff[ind]))




 
    seds_disk = seds_bands[0]
    seds_bulge = seds_bands[1]
    seds_total = seds_bands[2]

    L1p4radio = seds_total[12,:]

    #compute luminosity function
    ind = np.where((lir_total > 0) & (lir_total < 1e20))
    H, bins_edges = np.histogram(np.log10(lir_total[ind]),bins=np.append(mbins,mupp))
    LFs_dust[index,:] = LFs_dust[index,:] + H

    #median Lradio vs LIR
    ind = np.where((lir_total > 1e6) & (lir_total < 1e20))
    meds_radio = bin_it(x=lir_total[ind], y=L1p4radio[ind])
  
    return(volh, h0)
    
def plot_lir_lf(plt, outdir, obsdir, LFs_dust, file_name):

    fig = plt.figure(figsize=(6,6))
    ytit = "$\\rm log_{10} (\\rm \\phi/\\, cMpc^{-3}\\, dex^{-1})$"
    xtit = "$\\rm log_{10} (L_{\\rm TIR}/L_{\\odot})$"
    xmin, xmax, ymin, ymax = 8.5, 12.8, -6, -2
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    cols = ('Indigo','purple','Navy','MediumBlue','Green','MediumAquamarine','LightGreen','YellowGreen','Gold','Orange','Coral','OrangeRed','red','DarkRed','FireBrick','Crimson','IndianRed','LightCoral','Maroon','brown','Sienna','SaddleBrown','Chocolate','Peru','DarkGoldenrod','Goldenrod','SandyBrown')
    ax = fig.add_subplot(111)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.5, 0.5, 0.5, 0.5))

    for b in range(0,len(zlist)):
        inp = np.where(LFs_dust[b,:] != 0)
        x = xmf[inp]
        y = LFs_dust[b,inp]
        ax.plot(x, y[0], linestyle='solid',color=cols[b], label='z=%s' % str(zlist[b]))
        print('#redshift: %s' % str(zlist[b]))
        for a,b in zip(x,y[0]):
            print (a,b)
    common.prepare_legend(ax, cols, loc=3)
    common.savefig(outdir, fig, 'LIR_total_highz_'+file_name+'.pdf')



def main(model_dir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()

    file_name = "eagle-rr14-radio-only"
    file_hdf5_sed = "Shark-SED-" + file_name + ".hdf5"

    fields_sed = {'SED/lir_dust': ('disk','bulge_t','total'),}
    fields_seds_bc = {'SED/lir_dust_contribution_bc': ('disk','bulge_t','total'),}
    fields_seds_bands = {'SED/ab_dust':('disk','bulge_t','total'),}

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge','sfr_disk', 'sfr_burst')}

    #Bands information:
    #(0): "z_SDSS", "Band_ionising_photons", "FUV_Nathan", "Band9_ALMA",
    #(4): "Band8_ALMA", "Band7_ALMA", "Band6_ALMA", "Band4_ALMA", "Band3_ALMA",
    #(9): "BandX_VLA", "BandC_VLA", "BandS_VLA", "BandL_VLA", "Band_610MHz",
    #(14): "Band_325MHz", "Band_150MHz"
 
    LFs_dust     = np.zeros(shape = (len(zlist), len(mbins)))

    for index, snapshot in enumerate(redshift_table[zlist]):
        print("Will read snapshot %s" % (str(snapshot)))
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        seds_lir = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed)
        seds_bands = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_seds_bands, subvols, file_hdf5_sed)
        seds_lir_bc = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_seds_bc, subvols, file_hdf5_sed)


        (volh, h0) = prepare_data(hdf5_data, seds_lir, seds_bands, seds_lir_bc, index, LFs_dust, obsdir)

    def take_log(x,v,h):
        x = x / (v / h**3.0)
        ind = np.where(x > 0)
        x[ind] = np.log10(x[ind])
        return x

    LFs_dust = take_log(LFs_dust, volh, h0)
    plot_lir_lf(plt, outdir, obsdir, LFs_dust, file_name)

if __name__ == '__main__':
    main(*common.parse_args())
