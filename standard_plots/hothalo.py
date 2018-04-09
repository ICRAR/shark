#
#    ICRAR - International Centre for Radio Astronomy Research
#    (c) UWA - The University of Western Australia, 2018
#    Copyright by UWA (in the framework of the ICRAR)
#    All rights reserved
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#    MA 02111-1307  USA
#
"""Hot halo plots"""

import numpy as np

import common
import utilities_statistics as us

##################################
# Constants
ConvKtoeV = 1.0 / (1.160452e4)
lcool_sim_to_cgs = 6.3031831e-14

tlow = -2.0
tupp = 2.0
dm = 0.15
tbins = np.arange(tlow,tupp,dm)
tfunc = tbins + dm/2.0

def prepare_data(hdf5_data):

    typeg, vhalo, cooling_rate = hdf5_data
    tvir = np.zeros(shape = len(typeg))
    Lcool = np.zeros(shape = len(typeg))

    ind = np.where((typeg == 0) & (cooling_rate > 0) & (vhalo > 0))
    tvir[ind]  = 97.48 * pow(vhalo[ind],2.0) * ConvKtoeV / 1e3 #in keV
    Lcool[ind] = 0.5  * cooling_rate[ind] * pow(vhalo[ind] ,2.0) * lcool_sim_to_cgs #cooling luminosity in units of 10^40 erg/s

    return us.wmedians(x=np.log10(tvir[ind]), y=np.log10(Lcool[ind]), xbins=tfunc)

def plot_cooling_rate(plt, outdir, med_tvir):

    fig = plt.figure(figsize=(5,5))
    xmin, xmax, ymin, ymax = -1.8, 1.2, -1, 6

    ax = fig.add_subplot(111)
    xtit="$\\rm log_{10}(T_{\\rm vir}/\\rm keV)$"
    ytit="$\\rm log_{10}(L_{\\rm cool}/ 10^{40}\\rm erg\,s^{-1})$"
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    xleg= xmax - 0.2 * (xmax - xmin)
    yleg= ymax - 0.1 * (ymax - ymin)
    ax.text(xleg, yleg, 'z=0')

    ind = np.where(med_tvir[0, :] != 0)
    xplot = tfunc[ind]
    yplot = med_tvir[0, ind]
    errdn = med_tvir[1, ind]
    errup = med_tvir[2, ind]

    ax.errorbar(xplot,yplot[0],color='k', label="SHArk")
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='+',markersize=2)
    ax.plot(tfunc, 3.0 * tfunc + 1.9, 'r', linestyle='dashed', label='Anderson+15')

    common.prepare_legend(ax, ['r','k'], loc=2)
    common.savefig(outdir, fig, 'cooling_rate.pdf')

def main():

    plt = common.load_matplotlib()
    modeldir, outdir, snapshot = common.parse_args(requires_observations=False)

    fields = {'Galaxies': ('type', 'vvir_hosthalo', 'cooling_rate')}
    hdf5_data = common.read_data(modeldir, snapshot, fields, include_h0_volh=False)
    med_tvir = prepare_data(hdf5_data)

    plot_cooling_rate(plt, outdir, med_tvir)

if __name__ == '__main__':
    main()
