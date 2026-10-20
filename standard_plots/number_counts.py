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
"""Number counts"""

import functools

import numpy as np

import common
import utilities_statistics as us

from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
import astropy.units as u
import astropy.cosmology.units as cu


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

def prepare_data(seds_bands, index, nc_dust, obsdir):

    lir_total = seds_bands[1] #total apparent magnitude with dust

    #calculate luminosity distance for this redshift
    fluxes = 10**((seds_total[:,:] - 8.9) / (-2.5)) * 1e3 #in mJy

    #calculate the area of the survey at the given redshift
    dL = cosmology.funcs.angular_diameter_distance(zlist[index] * cu.redshift) / 57.2958 #Mpc per 1 degree
    area = ((volh / h0**3)**(1/3) / dL) **2 #in degree^2 

    #calculate the delta redshift at the given redshift
    dc = cosmo.comoving_distance(zlist[index] * cu.redshift)
    d1 = dc - 0.5 * (volh / h0**3)**(1/3) * u.Mpc
    d2 = dc + 0.5 * (volh / h0**3)**(1/3) * u.Mpc
    r1 = d1.to(cu.redshift, cu.redshift_distance(cosmo, kind="comoving", zmax=1200))
    r2 = d2.to(cu.redshift, cu.redshift_distance(cosmo, kind="comoving", zmax=1200))
    delta_z = r2 - r1 #given the length of the box, what's the delta_z.

    #and now we build the number counts per area per delta_redshift



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

    fields_seds_bands = {'SED/ab_dust':('disk','total'),}

    #Bands information:
    #(0): "z_SDSS", "Band_ionising_photons", "FUV_Nathan", "Band9_ALMA",
    #(4): "Band8_ALMA", "Band7_ALMA", "Band6_ALMA", "Band4_ALMA", "Band3_ALMA",
    #(9): "BandX_VLA", "BandC_VLA", "BandS_VLA", "BandL_VLA", "Band_610MHz",
    #(14): "Band_325MHz", "Band_150MHz"
 
    nc_dust_z     = np.zeros(shape = (len(zlist), len(fbins)))
    nc_dust       = np.zeros(shape = (len(fbins)))

    for index, snapshot in enumerate(redshift_table[zlist]):
        print("Will read snapshot %s" % (str(snapshot)))
        seds_bands = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_seds_bands, subvols, file_hdf5_sed)

        (volh, h0) = prepare_data(seds_bands, index, nc_dust_z, obsdir)

    plot_lir_lf(plt, outdir, obsdir, LFs_dust, file_name)

if __name__ == '__main__':
    main(*common.parse_args())
