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


zlist=np.array([0,0.1,0.2,0.25,0.3,0.35,0.4,0.5,0.75,1.0])

##################################
#Constants
RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun
c_light  = 299792458.0 #m/s
PI       = 3.141592654
zsun = 0.0189
MpctoKpc = 1e3

mlow = 9.5
mupp = 12.5
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

mrlow = -6.0
mrupp = -1.0
dmr = 0.2
mrbins = np.arange(mrlow,mrupp,dmr)
xmrf = mrbins + dmr/2.0



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


def prepare_data(hdf5_data, index, dust_mass_z, dust_ratio_dist_z):

    #read properties from hdf5 file
    (h0, volh, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     sfr_disk, sfr_bulge, typeg,  mgas_disk, mgas_bulge, matom_disk, mmol_disk, matom_bulge, mmol_bulge, mvir_hosthalo, idtree, mzd, mzb) = hdf5_data
   
    bin_it   = functools.partial(us.wmedians, xbins=xmf)
 
    mstar = np.zeros(shape = len(mdisk))
    ind = np.where(mdisk + mbulge > 0)
    mstar[ind] = np.log10((mdisk[ind] + mbulge[ind])/h0)

    thresh= 10.8
    (mdustd, DToM_MW) = dust_mass(mzd, mgas_disk, h0)
    (mdustb, DToM_MW) = dust_mass(mzb, mgas_bulge, h0)

    #select MAGPI primary galaxies
    ind = np.where(mstar >= thresh)
    mmid, mhigh = np.percentile(np.log10(mvir_hosthalo[ind]), (33.333, 66.666))
    #select all MAGPI ``well-resolved'' sample
    ind = np.where((mstar >= thresh) | ((mstar > 9.5) & (typeg > 0)))
    mstar_prim = mstar[ind]
    print(len(mstar_prim))
    mdust_prim = (mdustd[ind] + mdustb[ind])
    mhalo_prim = np.log10(mvir_hosthalo[ind])
 
    #now study properties of dust in three environment densities
    ind = np.where((mhalo_prim < mmid) & (mdust_prim > 0))
    print('low density environment:', len(mstar_prim[ind]))
    dust_mass_z[index,0,:] = bin_it(x = mstar_prim[ind], y = np.log10(mdust_prim[ind]) - mstar_prim[ind])
    H_, _ = np.histogram(np.log10(mdust_prim[ind]) - mstar_prim[ind],bins=np.append(mrbins,mrupp))
    dust_ratio_dist_z[index,0,:] = dust_ratio_dist_z[index,0,:] + H_

    ind = np.where((mhalo_prim >= mmid) & (mhalo_prim < mhigh) & (mdust_prim > 0))
    dust_mass_z[index,1,:] = bin_it(x = mstar_prim[ind], y = np.log10(mdust_prim[ind]) - mstar_prim[ind])    
    H_, _ = np.histogram(np.log10(mdust_prim[ind]) - mstar_prim[ind],bins=np.append(mrbins,mrupp))
    dust_ratio_dist_z[index,1,:] = dust_ratio_dist_z[index,1,:] + H_

    ind = np.where((mhalo_prim >= mhigh) & (mdust_prim > 0))
    dust_mass_z[index,2,:] = bin_it(x = mstar_prim[ind], y = np.log10(mdust_prim[ind]) - mstar_prim[ind])    
    H_, _ = np.histogram(np.log10(mdust_prim[ind]) - mstar_prim[ind],bins=np.append(mrbins,mrupp))
    dust_ratio_dist_z[index,2,:] = dust_ratio_dist_z[index,2,:] + H_

    return(volh, h0)
    
def plot_dust_to_star_relation(plt, outdir, obsdir, dust_mass_z):

    fig = plt.figure(figsize=(11.5,4.5))
    xtit = "$\\rm log_{10}(M_{\\star}/M_{\\odot})$"
    ytit = "$\\rm log_{10}(M_{\\rm dust}/M_{\\star})$"
    xmin, xmax, ymin, ymax = 9.5, 12, -6, -2
    xleg = xmax - 0.8 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = [131, 132, 133]
    bands = ['low density', 'intermediate density', 'high density']
    cols = ['DarkRed','Red','DarkOrange','Gold','Chocolate','YellowGreen','LimeGreen','DarkTurquoise','LightSteelBlue','Navy']

    for i in range(0,len(subplots)):
        ax = fig.add_subplot(subplots[i])
        if i == 0:
           ytitle = ytit
        else:
           ytitle = ' '
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle, locators=(0.1, 1, 1, 1))
        ax.text(xleg,yleg,bands[i],fontsize=12)
        for z in range(0,len(zlist)):
            inp = np.where(dust_mass_z[z,i,0,:] != 0)
            yplot = dust_mass_z[z,i,0,inp]
            print(yplot.shape, z)
            ax.plot(xmf[inp], yplot[0,:], linestyle='solid',color=cols[z], label='z=%s' % str(zlist[z]))
    plt.subplots_adjust(right=0.87)
    common.prepare_legend(ax, cols, bbox_to_anchor=(1, 0., 0.5, 0.5))
    common.savefig(outdir, fig, 'DustToStellarMass_vs_redshift_environment.pdf')


    #plot MAGPI redshift only
    fig = plt.figure(figsize=(4,4.5))
    xmin, xmax, ymin, ymax = 9.5, 12, -6, -2
    bands = ['low density', 'intermediate density', 'high density']
    cols = ['Navy', 'Gold', 'Red']

    z = 3 #0.25
    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 1, 1))
    plt.subplots_adjust(left=0.15)

    for d in range(0,len(bands)):
        inp = np.where(dust_mass_z[z,d,0,:] != 0)
        yplot = dust_mass_z[z,d,0,inp]
        yplot_err_dn = dust_mass_z[z,d,1,inp]
        yplot_err_up = dust_mass_z[z,d,2,inp]
        ax.fill_between(xmf[inp],yplot[0]-yplot_err_dn[0],yplot[0]+yplot_err_up[0], facecolor=cols[d], alpha=0.5, interpolate=True)
        ax.plot(xmf[inp], yplot[0], linestyle='solid',color=cols[d], label=bands[d])
        print('Environment', d)
        for a,b,c,d in zip(xmf[inp],yplot[0],yplot_err_dn[0],yplot_err_up[0]):
            print(a,b,c,d)
    common.prepare_legend(ax, cols, loc='upper left')
    common.savefig(outdir, fig, 'DustToStellarMass_z0p25_environment.pdf')



def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'sfr_disk', 'sfr_burst', 'type', 
                           'mgas_disk', 'mgas_bulge','matom_disk', 'mmol_disk', 
                           'matom_bulge', 'mmol_bulge', 'mvir_hosthalo', 'id_halo_tree',  'mgas_metals_disk',
                           'mgas_metals_bulge')}

    dust_mass_z = np.zeros(shape = (len(zlist), 3, 3, len(mbins)))
    dust_ratio_dist_z = np.zeros(shape = (len(zlist), 3, len(mrbins)))

    for index, snapshot in enumerate(redshift_table[zlist]):
        print("Will read snapshot %s" % (str(snapshot)))
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)

        (volh, h0) = prepare_data(hdf5_data, index, dust_mass_z, dust_ratio_dist_z)

    plot_dust_to_star_relation(plt, outdir, obsdir, dust_mass_z)

if __name__ == '__main__':
    main(*common.parse_args())
