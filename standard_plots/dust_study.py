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

mlow = 7.0
mupp = 12.5
dm = 1.0
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

mrlow = -6.0
mrupp = -1.0
dmr = 0.2
mrbins = np.arange(mrlow,mrupp,dmr)
xmrf = mrbins + dmr/2.0

slow = -1
supp = 3.0
ds = 1.0
sbins = np.arange(slow,supp,ds)
xsf = sbins + ds/2.0

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


def prepare_data(hdf5_data, seds, seds_read, index, dust_mass_z, dust_ratio_dist_z, dust_sm_scaling, dust_sfr_scaling, s6_sm_scaling, s6_sfr_scaling, mean_sm_mode, mean_sfr_mode):

    #read properties from hdf5 file
    (h0, volh, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     sfr_disk, sfr_bulge, typeg,  mgas_disk, mgas_bulge, matom_disk, mmol_disk, matom_bulge, mmol_bulge, mvir_hosthalo, idtree, mzd, mzb) = hdf5_data



    if(seds_read):
       seds_tot = seds[1]
       band6 = 31
       sband6 = 10**((seds_tot[band6, :] - 8.9) / (-2.5)) #in Jy
       nu_obs = 249.827 #GHz
       nu_rest = nu_obs * (1 + zlist[index])
       Td = 25.0
       Tcmb = 2.73
       beta = 1.8
       Ttrue = (Td**(4+beta) + Tcmb**(4+beta) * ((1+zlist[index])**(4+beta) - 1 ))**(1.0/(4+beta))
       Tobs = Ttrue / (1 + zlist[index])
       nu0_nurest_fac = (352.6 / nu_rest)**beta
       temp_fact = 6.62607015e-34 / 1.380649e-23 * (nu_obs * 1e9)
   
       Bb_T = 1474.5 * (nu_obs)**3 / (np.exp(temp_fact / Tobs) - 1) #in Jy str^-1
      
       mdustv2 = 5.03e-31 * sband6 * (dl[index] * 3.086e22)**2 * nu0_nurest_fac / (1 + zlist[index])**4 / Bb_T / 0.0431
       
 
    bin_it   = functools.partial(us.wmedians, xbins=xmf)
    bin_it_sfr   = functools.partial(us.wmedians, xbins=xsf)

    ind = np.where(mdisk + mbulge > 0)
    
    mstar = np.log10((mdisk[ind] + mbulge[ind])/h0)
    sfr = np.log10((sfr_disk[ind] + sfr_bulge[ind])/1e9/h0)
    sfr_mode = sfr_bulge[ind]/(sfr_disk[ind] + sfr_bulge[ind])
    typeg = typeg[ind]
    mvir_hosthalo = mvir_hosthalo[ind]
    sfr_disk = sfr_disk[ind] 
    sfr_bulge = sfr_bulge[ind]
    (mdustd, DToM_MW) = dust_mass(mzd[ind], mgas_disk[ind], h0)
    (mdustb, DToM_MW) = dust_mass(mzb[ind], mgas_bulge[ind], h0)
    mdust = mdustd + mdustb

    if(seds_read):
       ind = np.where((mdust > 1e5) & (mstar > 8))
       print(np.median(mdust[ind]/mdustv2[ind]), np.std(mdust[ind]/mdustv2[ind]))
     

    ind = np.where((sfr_disk + sfr_bulge) <= 0)
    sfr[ind] = -10

    if(zlist[index] <= 2):
       delta_ms = (sfr - mstar) - (-11 + 0.5 * zlist[index])
    else:
       delta_ms = (sfr - mstar) - (-10)


    thresh= 10.

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

    mdustin = mdust
    ind = np.where((mstar >= 7) & (mdustin > 0))
    dust_sm_scaling[index,:] = us.stacking(x = mstar[ind], y = mdustin[ind], xbins=xmf)    
    dust_sfr_scaling[index,:] = us.stacking(x = sfr[ind], y = mdustin[ind], xbins=xsf)
    mean_sm_mode[index,:] = bin_it(x= mstar[ind], y=sfr_mode[ind])
    mean_sfr_mode[index,:] = bin_it_sfr(x= sfr[ind], y=sfr_mode[ind])

    if(seds_read):
       s6_sm_scaling[index,:] = us.stacking(x = mstar[ind], y = sband6[ind], xbins=xmf)
       s6_sfr_scaling[index,:] = us.stacking(x = sfr[ind], y = sband6[ind], xbins=xsf)

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

def plot_dust_to_star_relation_zevol(plt, outdir, obsdir, dust_mass_z, dust_sfr_scaling, s6_sm_scaling, s6_sfr_scaling, mean_sm_mode, mean_sfr_mode):

    fig = plt.figure(figsize=(11.5,4.5))
    xtit = "$\\rm log_{10}(M_{\\star}/M_{\\odot})$"
    ytit = "$\\rm log_{10}(M_{\\rm dust}/M_{\\odot})$"
    xmin, xmax, ymin, ymax =8, 12, 5, 8.5
    xleg = xmax - 0.8 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = [121, 122]
    cols = ['DarkRed','Red','DarkOrange','Gold','Chocolate','YellowGreen','LimeGreen','DarkTurquoise','LightSteelBlue','Navy']

    zlist_indx = [3,4,5,6,7,8,9,10,11,12] #,13] #,14]
    ax = fig.add_subplot(subplots[0])
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 1, 1))
    for z in range(0,len(zlist_indx)):

        inp = np.where(dust_mass_z[zlist_indx[z],1,:] != 0)
        xplot = dust_mass_z[zlist_indx[z],0,inp]
        yplot = dust_mass_z[zlist_indx[z],1,inp]
        mode_med = mean_sm_mode[zlist_indx[z],0,inp]
        zplot = s6_sm_scaling[zlist_indx[z],1,inp]
        ax.plot(xplot[0,:], yplot[0,:], linestyle='solid',color=cols[z], label='z=%s' % str(zlist[zlist_indx[z]]))
        print("#Mean Mstar, Mdust, S6, SFR_burst/SFR_total at redshift:", zlist[zlist_indx[z]])
        for a,b,c,d in zip(xplot[0,:], yplot[0,:], zplot[0,:], mode_med[0,:]):
            print(a,b,c,d)

    plt.subplots_adjust(right=0.87)
    common.prepare_legend(ax, cols, loc=4)


    xtit = "$\\rm log_{10}(\\rm SFR/M_{\\odot}\\, yr^{-1})$"
    ytit = "$\\rm log_{10}(M_{\\rm dust}/M_{\\odot})$"
    xmin, xmax, ymin, ymax =-1.5, 3, 5, 8.5

    ax = fig.add_subplot(subplots[1])
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 1, 1))
    for z in range(0,len(zlist_indx)):
        inp = np.where(dust_sfr_scaling[zlist_indx[z],1,:] != 0)
        xplot = dust_sfr_scaling[zlist_indx[z],0,inp]
        yplot = dust_sfr_scaling[zlist_indx[z],1,inp]
        mode_med = mean_sfr_mode[zlist_indx[z],0,inp]
        zplot = s6_sfr_scaling[zlist_indx[z],1,inp]
        ax.plot(xplot[0,:], yplot[0,:], linestyle='solid',color=cols[z])
        print("#Mean SFR, Mdust, S6, SFR_burst/SFR_total at redshift:", zlist[zlist_indx[z]])
        for a,b,c,d in zip(xplot[0,:], yplot[0,:], zplot[0,:], mode_med[0,:]):
            print(a,b,c,d)

    common.savefig(outdir, fig, 'DustToStellarMass_Scaling_vs.pdf')

def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'sfr_disk', 'sfr_burst', 'type', 
                           'mgas_disk', 'mgas_bulge','matom_disk', 'mmol_disk', 
                           'matom_bulge', 'mmol_bulge', 'mvir_hosthalo', 'id_halo_tree',  'mgas_metals_disk',
                           'mgas_metals_bulge')}

    fields_sed = {'SED/ap_dust': ('disk','total'),}
    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5"

    dust_mass_z = np.zeros(shape = (len(zlist), 3, 3, len(mbins)))
    dust_ratio_dist_z = np.zeros(shape = (len(zlist), 3, len(mrbins)))
    dust_sm_scaling =  np.zeros(shape = (len(zlist), 2, len(mbins)))
    dust_sfr_scaling =  np.zeros(shape = (len(zlist), 2, len(sbins)))
    mean_sm_mode =  np.zeros(shape = (len(zlist), 3, len(mbins)))
    mean_sfr_mode =  np.zeros(shape = (len(zlist), 3, len(sbins)))
    s6_sm_scaling =  np.zeros(shape = (len(zlist), 2, len(mbins)))
    s6_sfr_scaling =  np.zeros(shape = (len(zlist), 2, len(sbins)))

    #(0): "FUV_GALEX", "NUV_GALEX", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS",
    #(6): "z_SDSS", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1_WISE",
    #(12): "I1_Spitzer", "I2_Spitzer", "W2_WISE", "I3_Spitzer", "I4_Spitzer",
    #(17): "W3_WISE", "W4_WISE", "P70_Herschel", "P100_Herschel",
    #(21): "P160_Herschel", "S250_Herschel", "S350_Herschel", "S450_JCMT",
    #(25): "S500_Herschel", "S850_JCMT", "FUV_Nathan", "Band9_ALMA",
    #(29): "Band8_ALMA", "Band7_ALMA", "Band6_ALMA", "Band4_ALMA"

    seds_read = True

    for index, snapshot in enumerate(redshift_table[zlist]):
        print("Will read snapshot %s" % (str(snapshot)), " corresponding to redshift,", zlist[index])
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        if(seds_read):
           seds = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed, subvols, file_hdf5_sed)
        else:
           seds = []

        (volh, h0) = prepare_data(hdf5_data, seds, seds_read, index, dust_mass_z, dust_ratio_dist_z, dust_sm_scaling, dust_sfr_scaling,
                                 s6_sm_scaling, s6_sfr_scaling, mean_sm_mode, mean_sfr_mode)

    #plot_dust_to_star_relation(plt, outdir, obsdir, dust_mass_z)
    plot_dust_to_star_relation_zevol(plt, outdir, obsdir, dust_sm_scaling, dust_sfr_scaling, s6_sm_scaling, s6_sfr_scaling, mean_sm_mode, mean_sfr_mode)

if __name__ == '__main__':
    main(*common.parse_args())
