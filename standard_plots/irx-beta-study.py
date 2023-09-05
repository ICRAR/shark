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


zlist=np.array([0.194739]) #, 0.9, 2.0, 3.0, 3.95, 5, 5.96, 7, 8, 8.94]) #, 8.94]) #, 9.95])

#0.194739, 0.254144, 0.359789, 0.450678, 0.8, 0.849027, 0.9, 1.20911, 1.28174, 1.39519, 1.59696, 2.00392, 2.47464723643932, 2.76734390952347, 3.01916, 3.21899984389701, 3.50099697082904, 3.7248038025221, 3.95972, 4.465197621546, 4.73693842543988, 5.02220991014863, 5.2202206934302, 5.52950356184419, 5.74417977285603, 5.96593, 6.19496927748119, 6.55269895697227, 7.05756323172746, 7.45816170313544, 7.73629493731708, 8.02352,8.32018565809831, 8.47220854014322, 8.78358705435761, 8.94312532315157, 9.27010372804765, 9.437541750167, 9.78074128377067, 9.95655])

##################################
#Constants
RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun
c_light  = 299792458.0 #m/s
PI       = 3.141592654

mlow = 6.5
mupp = 12.5
dm = 0.25
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

dmobs = 0.4
mbins_obs = np.arange(mlow,mupp,dmobs)
xmf_obs = mbins_obs + dmobs/2.0

vlow = -4.0
vupp = 3.0
dv   = 0.1
vbins = np.arange(vlow,vupp,dv)
xv    = vbins + dv/2.0

vlow2 = -4.0
vupp2 = 3.0
dv2   = 0.25
vbins2 = np.arange(vlow2,vupp2,dv2)
xv2    = vbins2 + dv2/2.0


btbins = [0, 0.05, 0.95, 1]

btlow = 0.0
btupp = 1.0
dbt   = 0.05
btbins2 = np.arange(btlow,btupp,dbt)
xbt    = btbins2 + dbt/2.0

zsun = 0.0189
#choose dust model between mm14, rr14 and constdust
m14 = False
rr14 = False
constdust = False
rr14xcoc = True


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

def smooth(x, y, ndeg):
    fit = np.polyfit(x,y,ndeg)
    print(fit) 
    y_smooth = np.zeros(shape = len(x))
    for j in range(0,ndeg+1):
        y_smooth[:] = y_smooth[:] + fit[j] * x[:]**(ndeg -j)

    return y_smooth

def prepare_data(hdf5_data, ext_data, seds, seds_lir, index, zsnap, obsdir, sm_irx):

    #read properties from hdf5 file
    (h0, volh, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     sfr_disk, sfr_bulge, typeg,  mgas_disk, mgas_bulge, matom_disk, mmol_disk, matom_bulge, mmol_bulge, mvir_hosthalo, 
     idtree, mzd, mzb) = hdf5_data

    (eta_bulge, eta_disk, tau_birth_bulge, tau_birth_disk, tau_screen_bulge, tau_screen_disk) = ext_data 

    sim_size = volh * h0**3.0 
    #compute dust masses
    (mdustd, DToM_MW) = dust_mass(mzd, mgas_disk, h0)
    (mdustb, DToM_MW) = dust_mass(mzb, mgas_bulge, h0)
    mdust = mdustb + mdustd

    #gas metallicity
    zgas = (mzd + mzb) / (mgas_disk + mgas_bulge) / zsun

    #star formation rate (in Msun/Gyr)
    sfr = (sfr_disk + sfr_bulge)/h0

    #define total stellar mass, bulde-to-total stellar mass ratio and read in SED files
    mstars_tot = (mdisk+mbulge)/h0
    bt = mbulge / (mdisk+mbulge)

    ind = np.where(mstars_tot > 0)
    zgas = zgas[ind]
    mstars = mstars_tot[ind]
    sfr = sfr[ind]
    typeg = typeg[ind]
    mdust= mdust[ind]
    SEDs_dust_total = seds[4] #total absolute magnitudes with dust
    lir_total = seds_lir[1] #total IR luminosity
    lir_disk = seds_lir[0] #disk IR luminosity
    lir_total = lir_total[0,:]
    lir_disk = lir_disk[0,:]

    #compute main sequence fit:
    ind = np.where((mstars > 1e9) & (mstars < 1e10) & (sfr > 0) & (typeg == 0))
    fit_ms = np.polyfit(np.log10(mstars[ind]), np.log10(sfr[ind]/mstars[ind]), 1)

    #select active galaxies
    active = np.zeros(shape = len(zgas))
    
    if zsnap <= 2:
       act = np.where(np.log10(sfr/mstars) > -2.0 + 0.5 * zsnap)
       active[act] = 1.0
    else:
       act = np.where(np.log10(sfr/mstars) > -1)
       active[act] = 1.0
   
    pos = np.where((SEDs_dust_total[0,:] > -30) & (SEDs_dust_total[0,:] < -10) & (mstars > 1e7) & (zgas > 0) & (typeg == 0))
    lfuv_total = 10.0**(SEDs_dust_total[0,pos]/(-2.5)) * 3631.0 * 0.06248261860502975 #in Lsun (for 1500 Angstrom)
    lnuv_total = 10.0**(SEDs_dust_total[1,pos]/(-2.5)) * 3631.0 * 0.03748957116301785 #in Lsun (for 2500 Angstrom)
    luv_rat = 10.0**(SEDs_dust_total[0,pos]/(-2.5)) / 10.0**(SEDs_dust_total[1,pos]/(-2.5)) #the ratio is on flux, not luminosity
    irx = lir_total[pos] / lfuv_total
    beta = np.log10(luv_rat) / (-0.22184874961635637) - 2.0 # number -0.22184874961635637 = np.log10(1500/2500)
    mstar_in = np.log10(mstars[pos])
    zgas_in = np.log10(zgas[pos])
    ssfr_in = np.log10(sfr[pos] / mstars[pos]) - (fit_ms[0] * np.log10(mstars[pos]) + fit_ms[1])
    sfr_in = np.log10(sfr[pos]/1e9/h0)
    mdust_in = np.log10(mdust[pos])
    eta_comp = (eta_bulge[pos] * (lir_total[pos] - lir_disk[pos]) +  eta_disk[pos] * lir_disk[pos]) / lir_total[pos] #luminosity weighted eta
    tau_birth_comp = (tau_birth_bulge[pos] * (lir_total[pos] - lir_disk[pos]) + tau_birth_disk[pos] * lir_disk[pos]) / lir_total[pos] #luminosity weighted tau_birth
    tau_screen_comp = (tau_screen_bulge[pos] * (lir_total[pos] - lir_disk[pos]) + tau_screen_disk[pos] * lir_disk[pos]) / lir_total[pos] #luminosity weighted tau_birth

    irx = irx[0,:]
    beta = beta[0,:] 
    pos = np.where(irx != 0)
    irx = np.log10(irx[pos])
    beta = beta[pos]
    mstar_in = mstar_in[pos]
    zgas_in = zgas_in[pos]
    ssfr_in = ssfr_in[pos]
    sfr_in = sfr_in[pos]
    mdust_in = mdust_in[pos]
    eta_comp = eta_comp[pos]
    tau_birth_comp = tau_birth_comp[pos]
    tau_screen_comp = tau_screen_comp[pos]

    sm_irx[index,:] = us.wmedians(x=mstar_in, y=irx, xbins=xmf)

    return(volh, h0, irx, beta, mstar_in, zgas_in, ssfr_in, mdust_in, eta_comp, tau_birth_comp, tau_screen_comp, sfr)
    
def plot_irx_beta(plt, outdir, obsdir, irx, beta, mstar, zgas_in, ssfr, mdust, eta, tau_birth, sfr, tau_screen, z, dust_model):

    ind = np.where((mstar >= 7.5) & (mstar <= 8.5))
    print("#Data for galaxies with Mstar/Msun in [7.5,8.5] for redshift:", z)
    print("#IRX beta log(Mstar/Msun) log(SFR/Msun yr^-1) Delta_MainSeq")
    for (a,b,c,d,e) in zip (irx[ind], beta[ind], mstar[ind], sfr[ind], ssfr[ind]):
        print(a,b,c,d,e)

    fig = plt.figure(figsize=(5,5))
    xtit = "$\\beta$"
    ytit = "log$_{10}(\\rm IRX)$"
    xmin, xmax, ymin, ymax = -2.5, 1.5, -1, 4.0
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15, bottom=0.17)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 0.5, 0.5))
    ax.text(xleg, yleg, 'z=%s' % str(z))
 
    ax.hexbin(beta, irx, xscale='linear', yscale='linear', gridsize=(50,50), cmap='jet', mincnt=10) 

    def plot_obs(ax):
        #plot attenuation curves
        #yv = np.log10((10.0**(0.4*(1.99*(vx + 2.23))) - 1.0) * 1.7)
        #ax.plot(xv, yv, linestyle='solid', color='k', label='Meurer+99')
        #consensus SMC from Bowens et al. (2020)
        yv = np.log10((10.0**(0.4*(1.1*(xv + 1.85))) - 1.0) * 1.7)
        ax.plot(xv, yv, linestyle='dotted', color='k', label='SMC')
        #Reddy+15
        A = 1.84 * xv + 4.48
        yv = np.log10(10.0**(0.4 * A) - 1.0) + 0.076
        ax.plot(xv, yv, linestyle='dashed', color='k', label='Reddy+15')
        #consensus z=0 from Bowens et al. (2020)
        yv = np.log10((10.0**(0.4*(1.86*(xv + 1.85))) - 1.0) * 1.7)
        ax.plot(xv, yv, linestyle='dashed', color='k', label='z=0')

    plot_obs(ax)
    common.prepare_legend(ax, 'k', loc=2)
    common.savefig(outdir, fig, 'irx_beta_%s_%s.pdf' % (str(z), dust_model))


    #plot IRX-beta for different stellar mass bins
    mstar_bins = [8,9,10,11,12]
    subplots = [141, 142, 143, 144]

    def plot_properties_in_irx_beta(prop, name, label, binsx=30, binsy=30, vmin=-1.5, vmax=0.5):
       fig = plt.figure(figsize=(11,3.7))
       xleg = xmin + 0.1 * (xmax - xmin)
       yleg = ymax - 0.35 * (ymax - ymin)
       xleg2 = xmin + 0.15 * (xmax - xmin)
       yleg2 = ymax + 0.05 * (ymax - ymin)
   
       meds_errors = np.zeros(shape = (4, 3, len(xv2)))
       for i,subplot in enumerate(subplots):
           ax = fig.add_subplot(subplot)
           if i == 0:
              ytitle = ytit
           else:
              ytitle = ''
           common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle, locators=(1, 1, 1, 1))
           if i == 0:
              ax.text(xleg, yleg, 'z=%s' % str(z))
           ax.text(xleg2, yleg2, '%s<$\\rm log_{10}(M_{\\star}/M_{\\odot})$<%s' % (str(mstar_bins[i]), str(mstar_bins[i+1])))
           inp = np.where((mstar > mstar_bins[i]) & (mstar <= mstar_bins[i+1]))
           x = beta[inp]
           y = irx[inp]
           meds = us.wmedians(x=x, y=y, xbins=xv2)
           meds_errors[i] = meds
           im = ax.hexbin(x, y, prop[inp], xscale='linear', yscale='linear', gridsize=(binsx,binsy), vmin=vmin, vmax=vmax, cmap='magma', mincnt=10, reduce_C_function=np.mean)
           pos = np.where(meds[0,:] != 0)
           y = meds[0,pos]
           yerrdn = meds[1,pos]
           yerrup = meds[2,pos]
           ax.errorbar(xv2[pos],y[0], yerr=[yerrdn[0],yerrup[0]], mfc='None', ecolor = 'DarkGray', mec='DarkGray', color='DarkGray', marker='o',linestyle='solid')
           #if(i == 0):
           #   print("#median IRX-beta relation for bin", mstar_bins[i], mstar_bins[i+1], " at redshift", z)
           #   for a, b, c, d in zip(xv2[pos],y[0], yerrdn[0],yerrup[0]):
           #       print(a,b,c,d)

           plot_obs(ax)
           if i == 0:
               common.prepare_legend(ax, 'k', loc=2)
           #else:
           #    plt.yticks([-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0],[" ", " "," "," "," "," "," "," "," "," "," "])
   
           if i == 3:
              cbar_ax = fig.add_axes([0.91, 0.17, 0.025, 0.7])
              cbar = fig.colorbar(im, cax=cbar_ax)
              cbar.ax.set_ylabel(label)
   
       plt.subplots_adjust(left=0.07, bottom=0.17)
       common.savefig(outdir, fig, 'irx_beta_stellar_mass_%s_%s_%s.pdf' % (name, str(z), dust_model))

       return meds_errors

    #plot metallicity, dust mass, and SSFR in irx-beta plane
    meds_errors = plot_properties_in_irx_beta(zgas_in, 'metallicity', '$\\rm log_{10}(Z_{\\rm gas}/Z_{\\odot})$', binsx=30, binsy=30, vmin=-1.5, vmax=0.5)
    meds_errors = plot_properties_in_irx_beta(ssfr, 'ssfr', '$\\rm log_{10}(sSFR/sSFR(MS))$', binsx=30, binsy=30, vmin=-1, vmax=1)
    meds_errors = plot_properties_in_irx_beta(mdust-mstar, 'mdust', '$\\rm log_{10}(M_{\\rm dust}/M_{\\star})$', binsx=30, binsy=30, vmin=-5, vmax=-2.5)

    #plot dust extinction parameters
    meds_errors = plot_properties_in_irx_beta(eta, 'eta', '$\\langle\\eta\\rangle$', binsx=30, binsy=30, vmin=-2.5, vmax=0)
    meds_errors = plot_properties_in_irx_beta(tau_birth, 'tau_birth', '$\\langle\\tau_{\\rm birth}\\rangle$', binsx=30, binsy=30, vmin=0, vmax=3)
    meds_errors = plot_properties_in_irx_beta(tau_screen, 'tau_screen', '$\\langle\\tau_{\\rm screen}\\rangle$', binsx=30, binsy=30, vmin=0, vmax=1.5)


    mstar_bins = [7,8,9,10,11,12.5]
    bins = 5
    #plot all stellar masses together
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15, bottom=0.17)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 0.5, 0.5))
    xleg = xmax - 0.25 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)
    ax.text(xleg, yleg, 'z=%s' % str(z))

    colors = ['CornflowerBlue', 'LightSeaGreen', 'Goldenrod', 'DarkOrange', 'DarkRed']
    meds_errors = np.zeros(shape = (5, 3, len(xv2)))

    for i in range(0,bins):
        inp = np.where((mstar > mstar_bins[i]) & (mstar <= mstar_bins[i+1]))
        x = beta[inp]
        y = irx[inp]
        meds = us.wmedians(x=x, y=y, xbins=xv2)
        meds_errors[i] = meds

        pos = np.where(meds_errors[i,0,:] != 0)
        y = meds_errors[i,0,pos]
        yerrdn = meds_errors[i,1,pos]
        yerrup = meds_errors[i,2,pos]
        ax.plot(xv2[pos],y[0], linestyle='solid', color=colors[i], label='%s<$\\rm log_{10}(M_{\\star}/M_{\\odot})$<%s' % (str(mstar_bins[i]), str(mstar_bins[i+1])))
        ax.fill_between(xv2[pos], y[0] - yerrdn[0], y[0] + yerrup[0], facecolor=colors[i], alpha=0.5, interpolate=True)
 
    plot_obs(ax)
    common.prepare_legend(ax, 'k', loc=2)
    common.savefig(outdir, fig, 'irx_beta_allmasses_%s_%s.pdf' % (str(z), dust_model))

def plot_irx_mass(plt, outdir, obsdir, sm_irx, zlist, dust_model):

    fig = plt.figure(figsize=(5,5))
    xtit = "log$_{10}(\\rm M_{\\star}/M_{\\odot})$"
    ytit = "log$_{10}(\\rm IRX)$"
    xmin, xmax, ymin, ymax = 7, 12, -2, 2.5
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    colors  = ('DarkRed','Salmon','Orange','YellowGreen','MediumSeaGreen','DarkTurquoise','MediumBlue','Purple','DarkSlateGray','Black')

    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15, bottom=0.17)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 0.5, 0.5))

    for i in range(0,len(zlist)):
        inp = np.where(sm_irx[i,0,:] != 0)
        y = sm_irx[i,0,inp]
        yerrdn = sm_irx[i,1,inp]
        yerrup = sm_irx[i,2,inp]
        ax.plot(xmf[inp], y[0,:], linestyle='solid', color=colors[i], label='z=%s' % str(zlist[i]))
        ax.fill_between(xmf[inp], y[0] - yerrdn[0], y[0] + yerrup[0], facecolor=colors[i], alpha=0.5, interpolate=True)

    common.prepare_legend(ax, colors, loc=4)
    common.savefig(outdir, fig, 'irx_stellarmass_allz_' + dust_model+ '.pdf')

def main(modeldir, outdir, redshift_table, subvols, obsdir):


    dust_model = "eagle-rr14"

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'sfr_disk', 'sfr_burst', 'type', 
                           'mgas_disk', 'mgas_bulge','matom_disk', 'mmol_disk', 
                           'matom_bulge', 'mmol_bulge', 'mvir_hosthalo', 'id_halo_tree', 'mgas_metals_disk', 'mgas_metals_bulge')}

    file_hdf5_sed = "Shark-SED-" + dust_model + ".hdf5"
    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total')}
    fields_sed_lir = {'SED/lir_dust': ('disk','total')}
    fields_ext = {'galaxies': ('pow_screen_bulge', 'pow_screen_disk', 'tau_birth_bulge', 'tau_birth_disk', 'tau_screen_bulge', 'tau_screen_disk')}

    #Bands information:
    #(0): "FUV_GALEX", "NUV_GALEX", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS",
    #(6): "z_SDSS", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1_WISE",
    #(12): "I1_Spitzer", "I2_Spitzer", "W2_WISE", "I3_Spitzer", "I4_Spitzer",
    #(17): "W3_WISE", "W4_WISE", "P70_Herschel", "P100_Herschel",
    #(21): "P160_Herschel", "S250_Herschel", "S350_Herschel", "S450_JCMT",
    #(25): "S500_Herschel", "S850_JCMT", "Band9_ALMA", "Band8_ALMA",
    #(29): "Band7_ALMA", "Band6_ALMA", "Band5_ALMA", "Band4_ALMA"

    #bands of interest band-7, band-6, band-4
    selec_alma = (29, 30, 32)
    flux_selec = (1e-10, 1e-2, 1e-1, 1.0) #to look at 0.01<S<0.1, 0.1<S<1, S>1mJy


    sm_irx = np.zeros(shape = (len(zlist), 3, len(xmf)))

    for index, snapshot in enumerate(redshift_table[zlist]):
        print("Will read snapshot %s" % (str(snapshot)))
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        seds = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed, subvols, file_hdf5_sed)
        seds_lir = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed_lir, subvols, file_hdf5_sed)
        ext_data = common.read_data_ext(modeldir, snapshot, fields_ext, subvols, dust_model)

        (volh, h0, irx, beta, mstar_in, zgas_in, ssfr, mdust, 
         eta, tau_birth, tau_screen, sfr) = prepare_data(hdf5_data, ext_data, seds, seds_lir, index, zlist[index], obsdir, sm_irx)
        plot_irx_beta(plt, outdir, obsdir, irx, beta, mstar_in, zgas_in, ssfr, mdust, eta, tau_birth, 
                      tau_screen, sfr, zlist[index], dust_model)

    def take_log(x,v,h):
        x = x / (v / h**3.0)
        ind = np.where(x > 0)
        x[ind] = np.log10(x[ind])
        return x

    plot_irx_mass(plt, outdir, obsdir, sm_irx, zlist, dust_model)

if __name__ == '__main__':
    main(*common.parse_args())
