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

#zlist=np.array([0,0.1,0.2,0.25,0.3,0.35,0.4,0.5,0.75,1.0, 1.5, 2.5, 3.5, 4.5, 5, 6])

zlist=np.array([0.194739, 0.450678, 0.8, 0.9, 1.20911, 1.59696, 2.00392, 2.47464723643932, 3.01916, 3.50099697082904, 3.95972, 4.465197621546, 5.02220991014863])#, 5.52950356184419]) #, 5.96593])
dl = np.array([981.08, 2576.79, 5159.47, 5962.96, 8582.07,  12092, 15970.7,  20643.1,  26237.9, 31322.9, 36259.8, 41791.3,  47981.7]) #,  53694.5]) #, 58659.5]) #in Mpc


##################################
#Constants

RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun
c_light  = 299792458.0 #m/s
PI       = 3.141592654
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


def prepare_data(hdf5_data, index, dust_density, redshifts):

    #read properties from hdf5 file
    (h0, volh, mgas_disk, mgas_bulge, mzd, mzb) = hdf5_data

    (mdustd, DToM_MW) = dust_mass(mzd, mgas_disk, h0)
    (mdustb, DToM_MW) = dust_mass(mzb, mgas_bulge, h0)
    mdust = mdustd + mdustb
    dust_density[index] = np.sum(mdust) / (volh / h0**3)
    print(redshifts[index], dust_density[index])

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mgas_disk', 'mgas_bulge', 'mgas_metals_disk',
                           'mgas_metals_bulge')}

    snap_data = np.loadtxt("/scratch/pawsey0119/clagos/medi-SURFS/1536/halos/trees/dhalo/hbt_trees_199/redshift_L210_N1536_list.txt")
    snapshots = range(40, 200, 1)
    snapshots_v2 = range(41, 201, 1)

    redshifts = snap_data[snapshots_v2[:],1]
    dust_density = np.zeros(shape = len(snapshots))


    for index, snapshot in enumerate(snapshots):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)

        prepare_data(hdf5_data, index, dust_density, redshifts)

if __name__ == '__main__':
    main(*common.parse_args())
