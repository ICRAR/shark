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


zlist=np.array([5.02220991014863, 5.52950356184419, 5.96593, 6.55269895697227, 7.05756323172746, 7.45816170313544, 8.02352, 9.95655])

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


def prepare_data(hdf5_data, seds, index, LFs_dust, obsdir):

    (h0, volh, mdisk, mbulge, sfr_disk, sfr_burst) = hdf5_data

    lir_total = seds[1] #total absolute magnitudes with dust

    #compute luminosity function
    ind = np.where((lir_total > 0) & (lir_total < 1e20))
    H, bins_edges = np.histogram(np.log10(lir_total[ind]),bins=np.append(mbins,mupp))
    LFs_dust[index,:] = LFs_dust[index,:] + H
  
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

    file_name = "eagle-rr14"
    file_hdf5_sed = "Shark-SED-" + file_name + ".hdf5"
    fields_sed = {'SED/lir_dust': ('disk','total'),}

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge','sfr_disk', 'sfr_burst')}

    #Bands information:
    #(0): "FUV_GALEX", "NUV_GALEX", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS",
    #(6): "z_SDSS", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1_WISE",
    #(12): "I1_Spitzer", "I2_Spitzer", "W2_WISE", "I3_Spitzer", "I4_Spitzer",
    #(17): "W3_WISE", "W4_WISE", "P70_Herschel", "P100_Herschel",
    #(21): "P160_Herschel", "S250_Herschel", "S350_Herschel", "S450_JCMT",
    #(25): "S500_Herschel", "S850_JCMT", "Band9_ALMA", "Band8_ALMA",
    #(29): "Band7_ALMA", "Band6_ALMA", "Band5_ALMA", "Band4_ALMA"

    LFs_dust     = np.zeros(shape = (len(zlist), len(mbins)))

    for index, snapshot in enumerate(redshift_table[zlist]):
        print("Will read snapshot %s" % (str(snapshot)))
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        seds = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed)
        (volh, h0) = prepare_data(hdf5_data, seds, index, LFs_dust, obsdir)

    def take_log(x,v,h):
        x = x / (v / h**3.0)
        ind = np.where(x > 0)
        x[ind] = np.log10(x[ind])
        return x

    LFs_dust = take_log(LFs_dust, volh, h0)
    plot_lir_lf(plt, outdir, obsdir, LFs_dust, file_name)

if __name__ == '__main__':
    main(*common.parse_args())
