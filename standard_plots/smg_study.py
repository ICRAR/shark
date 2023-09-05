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


zlist=np.array([0.194739, 0.254144, 0.359789, 0.450678, 0.8, 0.849027, 0.9, 1.20911, 1.28174, 1.39519, 1.59696, 2.00392, 2.47464723643932, 2.76734390952347, 3.01916, 3.21899984389701, 3.50099697082904, 3.7248038025221, 3.95972, 4.465197621546, 4.73693842543988, 5.02220991014863, 5.2202206934302, 5.52950356184419, 5.74417977285603, 5.96593, 6.19496927748119, 6.55269895697227, 7.05756323172746, 7.45816170313544, 7.73629493731708, 8.02352,8.32018565809831, 8.47220854014322, 8.78358705435761, 8.94312532315157, 9.27010372804765, 9.437541750167, 9.78074128377067, 9.95655])

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

def prepare_data(hdf5_data, seds, seds_nod, seds_ap, index, sfr_z, mmol_z, mdust_z, sfr_tot, mmol_tot, mdust_tot, flux_selec, selec_alma, frac_mvir_occupation, mvir_occupation, zsnap, obsdir, mstar_mf):

    #read properties from hdf5 file
    (h0, volh, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     sfr_disk, sfr_bulge, typeg,  mgas_disk, mgas_bulge, matom_disk, mmol_disk, matom_bulge, mmol_bulge, mvir_hosthalo, 
     idtree, mzd, mzb) = hdf5_data

    sim_size = volh / h0**3.0 
    #compute dust masses
    (mdustd, DToM_MW) = dust_mass(mzd, mgas_disk, h0)
    (mdustb, DToM_MW) = dust_mass(mzb, mgas_bulge, h0)
   
    #compute the total SFR,  molecular gas and dust masses in the box
    sfr_tot[index] = np.sum((sfr_disk + sfr_bulge) / 1e9 / h0)
    mmol_tot[index] = np.sum((mmol_bulge + mmol_disk) / h0)
    mdust_tot[index] = np.sum(mdustd + mdustb)

    #define total stellar mass, bulde-to-total stellar mass ratio and read in SED files
    mstars_tot = (mdisk+mbulge)/h0
    bt = mbulge / (mdisk+mbulge)
    ind = np.where(mstars_tot > 0)
    mstars = mstars_tot[ind]
    SEDs_dust_total = seds[4] #total absolute magnitudes with dust
    SEDs_vodust_total = seds_nod[1] #total absolute magnitudes no dust
    SEDs_app = seds_ap[4] #apparent magnitudes with dust
    A_nuv = SEDs_dust_total[1,:] - SEDs_vodust_total[1,:]
    negav = np.where(A_nuv < 0)
    A_nuv[negav] = 0.0

    color1 = SEDs_dust_total[1] - SEDs_dust_total[4]
    color2 = SEDs_dust_total[4] - SEDs_dust_total[8]

    quench_flag = np.zeros(shape = len(color1))
    quench = np.where((color1 > 3.0 * color2 + 1) & (color1 > 3.1))
    quench_flag[quench] = 1.0
    H, _ = np.histogram(np.log10(mstars[quench]), bins=np.append(mbins,mupp))
    mstar_mf[0,:,index+1] = mstar_mf[0,:,index+1] + H
    mstar_mf[0,:,index+1] = mstar_mf[0,:,index+1] / float(sim_size) / dm

    #select non-quenched galaxies
    non_quench = np.where(quench_flag == 0)
    H, _ = np.histogram(np.log10(mstars[non_quench]), bins=np.append(mbins,mupp))
    mstar_mf[1,:,index+1] = mstar_mf[1,:,index+1] + H
    mstar_mf[1,:,index+1] = mstar_mf[1,:,index+1] / float(sim_size) / dm

    # all galaxies
    H, _ = np.histogram(np.log10(mstars), bins=np.append(mbins,mupp))
    mstar_mf[2,:,index+1] = mstar_mf[2,:,index+1] + H
    mstar_mf[2,:,index+1] = mstar_mf[2,:,index+1] / float(sim_size) / dm

    #print some various properties of H-dropout galaxies (band = 9 is H-band and 13 is IRAC4.5microns)
    ind = np.where((SEDs_app[9,:] > 27) & (SEDs_app[13,:] < 24))
    print("Number of H-dropouts", len(mdisk[ind]))
    print("Volume density of H-dropouts", len(mdisk[ind])/(volh/h0**3.0))
    print("median stellar mass", np.median(mdisk[ind] + mbulge[ind])/h0)
    print("median SFR", np.median(sfr_disk[ind] + sfr_bulge[ind])/h0/1e9)

    #select galaxies with stellar masses > 0
    ind = np.where(mstars_tot > 0)
    sfrs_gals = (sfr_disk[ind] + sfr_bulge[ind]) / 1e9 / h0
    mmol_gals = (mmol_bulge[ind] + mmol_disk[ind])/h0
    mdust_gals = (mdustb[ind] + mdustd[ind])/h0
    types = typeg[ind]
    mvir = mvir_hosthalo[ind]
    idtrees = idtree[ind]

    #calculate total SFR of galaxies selected in different ALMA bands and different fluxes
    def calculate_lir(SEDs_dust_total, obsdir):
        file = obsdir+'/Models/Shark_SED_bands.dat'
        lambda_bands = np.loadtxt(file,usecols=[0],unpack=True)
        freq_bands   = c_light / (lambda_bands * 1e-10) #in Hz

        bands = [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 30] #from 8mu to 1000mu
        lir = np.zeros(shape = (len(SEDs_app[0,:])))
        fluxes = 10.0**(SEDs_dust_total[:,0,:] / -2.5) * 3631.0 * 1e-23 * 4.0 * PI * (3.086e19)**2.0 / 3.846e33 #in Lsun Hz−1
        print(fluxes.shape) 
        for i in range(0,len(fluxes[0,:])):
            for b in range(0,len(bands)-1):
                delta_freq = abs(freq_bands[bands[b]] - freq_bands[bands[b+1]])
                lir[i] = lir[i] + fluxes[bands[b], i] * delta_freq + abs(fluxes[bands[b+1], i] - fluxes[bands[b], i]) * delta_freq
        return lir
    ind = np.where(sfrs_gals > 50)
    lir_50sfr = calculate_lir(SEDs_dust_total[:,ind], obsdir)

    #compute the number density of bright IR galaxies (LIR > 3e12Lsun)
    ind = np.where(lir_50sfr > 3e12)
    numberdensity = (len(lir_50sfr[ind]) + 0.0) / (volh * h0**3.0)
    print('number density of galaxies with LIR>3e12Msun is', numberdensity, ' at redshift ', zsnap)

    #compute the total SFR and molecular gas mass contriubution of galaxies selected in different ALMA bands
    for b in range(0,len(selec_alma)):
        flux_gals_band = 10.0**(SEDs_app[selec_alma[b],:] / -2.5) * 3631.0 * 1e3 #in mJy
        for f in range(0,len(flux_selec)):
            if (f < 3):
                ind = np.where((flux_gals_band > flux_selec[f]) & (flux_gals_band <= flux_selec[f+1]))
            else:
                ind = np.where((flux_gals_band > flux_selec[f]) & (flux_gals_band < 1e10))
            sfr_z[b,f,index] = np.sum(sfrs_gals[ind])
            mmol_z[b,f,index] = np.sum(mmol_gals[ind])
            mdust_z[b,f,index] = np.sum(mdust_gals[ind])

    Calculate_ocuppation = False
    if(Calculate_ocuppation == True):
       print("number of unique halos", len(np.unique(idtrees)))
       #select most massive halos
       centrals = np.where(types == 0)
       print("number of halos with types==0", len(mvir[centrals]))
       mvir_centrals = mvir[centrals]
       id_mhalos = np.argsort(1.0/mvir_centrals) #sort from most massive to least massive
       ids_centrals = idtrees[centrals]
       idtree_sorted = ids_centrals[id_mhalos]
       mvir_sorted = mvir_centrals[id_mhalos]
       ids_mostmassive = idtree_sorted[0:20]
       mvir_occupation[index] = mvir_sorted[19]
       ngals = np.zeros(shape = (len(selec_alma)))
       for i in range(0,len(ids_mostmassive)):
           selec_group = np.where(idtrees == ids_mostmassive[i])
           if(len(idtrees[selec_group]) > 0):
              for b in range(0,len(selec_alma)):
                  flux_gals_band = 10.0**(SEDs_app[selec_alma[b],selec_group] / -2.5) * 3631.0 * 1e3 #in mJy
                  bright = np.where((flux_gals_band > max(flux_selec)) & (flux_gals_band < 1e10))
                  if(len(flux_gals_band[bright]) > 0):
                     ngals[b] = ngals[b] + 1
              
       for b in range(0,len(selec_alma)):
           frac_mvir_occupation[b,index] = ngals[b] / (len(ids_mostmassive) + 0.0)

    return(volh, h0)
    
def plot_sfr_contribution(plt, outdir, obsdir, sfr_z, sfr_tot, mmol_z, mmol_tot, mdust_z, mdust_tot, h0):

    #plot cosmic evolution of the SFR
    fig = plt.figure(figsize=(12,4.5))
    ytit = "$\\rm log_{10} (\\rm \\rho_{\\rm SFR}/ M_{\\odot} yr^{-1} cMpc^{-3})$"
    xtit = "redshift"
    xmin, xmax, ymin, ymax = 0, 10, -6, -1
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = [131, 132, 133]
    bands = ['band-7', 'band-6', 'band-4']
    cols = ['DarkBlue','MediumTurquoise','YellowGreen', 'Crimson']
    labels = ['$\\rm <0.01\\rm mJy$', '$\\rm  0.01-0.1\\rm mJy$', '$\\rm  0.1-1\\rm mJy$', '$\\rm  >1\\rm mJy$']

    def load_observations(ax, obsdir, h0):
        #Driver (Chabrier IMF), ['Baldry+2012, z<0.06']
        redD17d, redD17u, sfrD17, err1, err2, err3, err4 = common.load_observation(obsdir, 'Global/Driver18_sfr.dat', [0,1,2,3,4,5,6])
        hobs = 0.7
        xobsD17 = (redD17d+redD17u)/2.0
        yobsD17 = sfrD17 + np.log10(hobs/h0)
        errD17 = yobsD17*0. - 999.
        errD17 = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))
        ax.errorbar(xobsD17, yobsD17, yerr=[errD17,errD17], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o', label="Driver+18")

        redB12, sfrB12, errB12 = common.load_observation(obsdir, 'Global/Bouwens2012.dat', [0,1,2])
        hobs = 0.7
        yobsB12 = sfrB12 + np.log10(hobs/h0)
        ax.errorbar(redB12, yobsB12, yerr=[errB12,errB12], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s', label="Bouwens+12")

    for b in range(0,len(bands)):
        ax = fig.add_subplot(subplots[b])
        ytitle = ytit
        if(b > 0):
           ytitle=" "
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle, locators=(0.1, 1, 0.1, 1))
        ax.text(xleg,yleg,bands[b],fontsize=12) 
        for i in range(0,len(labels)):
            inp = np.where(sfr_z[b,i,:] != 0)
            x = zlist[inp]
            y = sfr_z[b,i,inp]
            print("SFR - Will print fitting parameters for band %s and selection %s" % (str(b), labels[i]))
            y_smooth = smooth(x, y[0], 3)
            ax.plot(x, y[0], linestyle='dotted',color=cols[i])
            ax.plot(x, y_smooth, linestyle='solid',color=cols[i], label=labels[i])
        ax.plot(zlist, sfr_tot,  linestyle='solid',color='k')
        load_observations(ax, obsdir, h0)
        if(b == 0):
           common.prepare_legend(ax, cols, loc=3)
    common.savefig(outdir, fig, 'SFR_evolution_SMG_contribution.pdf')


    fig = plt.figure(figsize=(12,2.5))
    ytit = "$\\rm log_{10}(fraction)$"
    xtit = "redshift"
    xmin, xmax, ymin, ymax = 0, 10, -2.5, 0
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymin + 0.1 * (ymax - ymin)

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
            y_smooth = smooth(zlist[inp], y[0], 3)
            y1 = sfr_z[b,i,inp] - sfr_tot[inp]
            y2 = y_smooth - sfr_tot[inp]
            ax.plot(zlist[inp], y1[0], linestyle='dotted',color=cols[i])
            ax.plot(zlist[inp], y2, linestyle='solid',color=cols[i])
        x=[0,10]
        y=[-1,-1]
        ax.plot(x,y,linestyle='dotted',color='k')
        fig.subplots_adjust(bottom=0.25)

    common.savefig(outdir, fig, 'fractional_SFR_evolution_SMG_contribution.pdf')

    # plot cosmic evolution of molecular gas
    fig = plt.figure(figsize=(12,4.5))
    ytit = "$\\rm log_{10} (\\rm \\rho_{\\rm H_2}/ M_{\\odot} cMpc^{-3})$"
    xtit = "redshift"
    xmin, xmax, ymin, ymax = -0.1, 10, 3.3, 8.3
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    def load_observations_h2(ax, obsdir, h0, caption=False):
        #Walter ASPECS ALMA program
        zloD16, zupD16, rhoH2loD16, rhoH2upD16  = common.load_observation(obsdir, 'Global/Decarli19_H2.dat', [0,1,2,3])
        zD16 =(zupD16 + zloD16)/2.0
        rhoH2D16 = (rhoH2loD16 + rhoH2upD16)/2.0
        hobs = 0.7
        xobs    = zD16
        errxlow = zD16-zloD16
        errxup  = zupD16-zD16
        yobs = rhoH2D16 + np.log10(pow(hobs/h0,3.0))
        errylow = rhoH2D16 - rhoH2loD16
        erryup  = rhoH2upD16 - rhoH2D16
        ax.errorbar(xobs, yobs, xerr=[errxlow,errxup], yerr=[errylow,erryup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='d',label="Decarli+19" if caption == True else None)

        #COLDz
        zloD16, zupD16, rhoH2loD16, rhoH2D16, rhoH2upD16  = common.load_observation(obsdir, 'Global/Riechers19_H2.dat', [0,1,2,3,4])
        zD16 =(zupD16 + zloD16)/2.0
        hobs = 0.7
        xobs    = zD16
        errxlow = zD16-zloD16
        errxup  = zupD16-zD16
        yobs = np.log10(rhoH2D16) + np.log10(pow(hobs/h0,3.0))
        errylow = np.log10(rhoH2D16) - np.log10(rhoH2loD16)
        erryup  = np.log10(rhoH2upD16) - np.log10(rhoH2D16)
        ax.errorbar(xobs, yobs, xerr=[errxlow,errxup], yerr=[errylow,erryup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Riechers+19" if caption == True else None)

        #ALMACAL-CO
        zloD16, zupD16, rhoH2loD16, rhoH2upD16  = common.load_observation(obsdir, 'Global/Hamanowicz20_H2.dat', [0,1,3,4])
        zD16 =(zupD16 + zloD16)/2.0
        rhoH2D16 = (rhoH2loD16 + rhoH2upD16)/2.0
        hobs = 0.7
        xobs    = zD16
        errxlow = zD16-zloD16
        errxup  = zupD16-zD16
        yobs = rhoH2D16 + np.log10(pow(hobs/h0,3.0))
        errylow = rhoH2D16 - rhoH2loD16
        erryup  = rhoH2upD16 - rhoH2D16
        ax.errorbar(xobs, yobs, xerr=[errxlow,errxup], yerr=[errylow,erryup], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='D',label="Hamanowicz+20" if caption == True else None)


        #z0 data
        zD16, zloD16, zupD16, rhoH2D16, rhoH2loD16, rhoH2upD16  = common.load_observation(obsdir, 'Global/H2_z0.dat', [0,1,2,3,4,5])
        xobs    = zD16
        errxlow = zD16-zloD16
        errxup  = zupD16-zD16
        yobs = np.log10(rhoH2D16) + np.log10(pow(hobs/h0,3.0))
        errylow = np.log10(rhoH2D16) - np.log10(rhoH2loD16)
        erryup  = np.log10(rhoH2upD16) - np.log10(rhoH2D16)
        ax.errorbar(xobs[0:1], yobs[0:1], xerr=[errxlow[0:1],errxup[0:1]], yerr=[errylow[0:1],erryup[0:1]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Boselli+14" if caption == True else None)
        ax.errorbar(xobs[1:2], yobs[1:2], xerr=[errxlow[1:2],errxup[1:2]], yerr=[errylow[1:2],erryup[1:2]], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*',label="Fletcher+20" if caption == True else None)

    for b in range(0,len(bands)):
        ax = fig.add_subplot(subplots[b])
        ytitle = ytit
        if(b > 0):
           ytitle=" "
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle, locators=(0.1, 1, 0.1, 1))
        ax.text(xleg,yleg,bands[b],fontsize=12) 
        for i in range(0,len(labels)):
            inp = np.where(mmol_z[b,i,:] != 0)
            x = zlist[inp]
            y = mmol_z[b,i,inp]
            print("H2 - Will print fitting parameters for band %s and selection %s" % (str(b), labels[i]))
            y_smooth = smooth(x, y[0], 3)
            ax.plot(x, y[0], linestyle='dotted',color=cols[i])
            ax.plot(x, y_smooth, linestyle='solid',color=cols[i], label=labels[i] if b == 0 else None)

        ax.plot(zlist, mmol_tot,  linestyle='solid',color='k')
        if(b == 1):
           load_observations_h2(ax, obsdir, h0, caption=True)
        else:
           load_observations_h2(ax, obsdir, h0, caption=False)
        if(b == 0):
           common.prepare_legend(ax, cols, loc=3)
        if(b == 1):
           common.prepare_legend(ax, ['k','k','k','k'], loc=3)

    common.savefig(outdir, fig, 'H2_evolution_SMG_contribution.pdf')


    fig = plt.figure(figsize=(12,2.5))
    ytit = "$\\rm log_{10}(fraction)$"
    xtit = "redshift"
    xmin, xmax, ymin, ymax = 0, 10, -2.5, 0
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymin + 0.1 * (ymax - ymin)

    for b in range(0,len(bands)):
        ax = fig.add_subplot(subplots[b])
        ytitle = ytit
        if(b > 0):
           ytitle=" "
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle, locators=(0.1, 1, 0.1, 1))
        ax.text(xleg,yleg,bands[b],fontsize=12) 
        for i in range(0,len(labels)):
            inp = np.where(mmol_z[b,i,:] != 0)
            y = mmol_z[b,i,inp]
            y_smooth = smooth(zlist[inp], y[0], 3)
            y1 = mmol_z[b,i,inp] - mmol_tot[inp]
            y2 = y_smooth - mmol_tot[inp]
            ax.plot(zlist[inp], y1[0], linestyle='dotted',color=cols[i])
            ax.plot(zlist[inp], y2, linestyle='solid',color=cols[i])
        x=[0,10]
        y=[-1,-1]
        ax.plot(x,y,linestyle='dotted',color='k')
        fig.subplots_adjust(bottom=0.25)

    common.savefig(outdir, fig, 'fractional_H2_evolution_SMG_contribution.pdf')

    #plot cosmic evolution of dust mass
    fig = plt.figure(figsize=(12,4.5))
    ytit = "$\\rm log_{10} (\\rm \\rho_{\\rm dust}/ M_{\\odot} cMpc^{-3})$"
    xtit = "redshift"
    xmin, xmax, ymin, ymax = -0.1, 5, 3., 6
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    for b in range(0,len(bands)):
        ax = fig.add_subplot(subplots[b])
        ytitle = ytit
        if(b > 0):
           ytitle=" "
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytitle, locators=(0.1, 1, 0.1, 1))
        ax.text(xleg,yleg,bands[b],fontsize=12) 
        for i in range(0,len(labels)):
            inp = np.where(mdust_z[b,i,:] != 0)
            x = zlist[inp]
            y = mdust_z[b,i,inp]
            y_smooth = smooth(x, y[0], 3)
            ax.plot(x, y[0], linestyle='dotted',color=cols[i])
            ax.plot(x, y_smooth, linestyle='solid',color=cols[i], label=labels[i] if b == 0 else None)

        ax.plot(zlist, mdust_tot,  linestyle='solid',color='k')
        print("Cosmic dust mass density")
        for a,b in zip(zlist, mdust_tot):
            print(a,b)
        if(b == 0):
           common.prepare_legend(ax, cols, loc=3)

    common.savefig(outdir, fig, 'dust_evolution_SMG_contribution.pdf')


def plot_occupation_massive_halos(plt, outdir, obsdir, frac_mvir_occupation):

    fig = plt.figure(figsize=(5,4.5))
    ytit = "$\\rm f_{\\rm occupation}$"
    xtit = "redshift"
    xmin, xmax, ymin, ymax = 0, 10, 0, 1
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = [131, 132, 133]
    bands = ['band-7', 'band-6', 'band-4']
    cols = ['MediumTurquoise','YellowGreen', 'Crimson']

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 0.1))
    ax.text(xleg,yleg,'S>1 mJy',fontsize=12) 

    for b in range(0,len(bands)):
        ax.plot(zlist, frac_mvir_occupation[b,:], linestyle='solid',color=cols[b], label=bands[b])
    common.prepare_legend(ax, cols, loc=3)
    common.savefig(outdir, fig, 'MassiveHalosContribution.pdf')

def plot_smf(plt, outdir, obsdir, mstar_mf):

    fig = plt.figure(figsize=(5,12))
    xtit = "log$_{10}$ M$_{*}$ [M$_{\odot}$]"
    ytit = "log$_{10}$ dn/dlog$_{10}$(M$_{*}$) [cMpc$^{-3}$]"
    xmin, xmax, ymin, ymax = 8, 12.5, -7, 0
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = [311, 312, 313]

    for s, subplot in enumerate(subplots):
        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.5, 0.5, 0.5, 0.5))
    
        for b in range(0,len(zlist)): 
            ax.plot(xmf, np.log10(mstar_mf[s,:,b+1]), linestyle='solid')
    #common.prepare_legend(ax, cols, loc=3)
    common.savefig(outdir, fig, 'smf_quiescent_sf.pdf')



def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'sfr_disk', 'sfr_burst', 'type', 
                           'mgas_disk', 'mgas_bulge','matom_disk', 'mmol_disk', 
                           'matom_bulge', 'mmol_bulge', 'mvir_hosthalo', 'id_halo_tree', 'mgas_metals_disk', 'mgas_metals_bulge')}

    file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5"
    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}
    fields_sed_nod = {'SED/ab_nodust': ('disk','total')}
    fields_sed_ap = {'SED/ap_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}

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

    mstar_mf = np.zeros(shape = (3, len(xmf), len(zlist) + 1))
    sfr_z = np.zeros(shape = (len(selec_alma), len(flux_selec), len(zlist)))
    mmol_z = np.zeros(shape = (len(selec_alma), len(flux_selec), len(zlist)))
    mdust_z = np.zeros(shape = (len(selec_alma), len(flux_selec), len(zlist)))
    sfr_tot = np.zeros(shape = (len(zlist)))
    mmol_tot = np.zeros(shape = (len(zlist)))
    mdust_tot = np.zeros(shape = (len(zlist)))
    frac_mvir_occupation = np.zeros(shape = (len(selec_alma), len(zlist)))
    mvir_occupation = np.zeros(shape = (len(zlist)))

    for index, snapshot in enumerate(redshift_table[zlist]):
        print("Will read snapshot %s" % (str(snapshot)))
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        seds = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed, subvols, file_hdf5_sed)
        seds_nod = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed_nod, subvols, file_hdf5_sed)
        seds_ap = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed_ap, subvols, file_hdf5_sed)

        (volh, h0) = prepare_data(hdf5_data, seds, seds_nod, seds_ap, index, sfr_z, mmol_z, mdust_z, sfr_tot, mmol_tot, mdust_tot, flux_selec, selec_alma, frac_mvir_occupation, mvir_occupation, zlist[index], obsdir, mstar_mf)

    def take_log(x,v,h):
        x = x / (v / h**3.0)
        ind = np.where(x > 0)
        x[ind] = np.log10(x[ind])
        return x

    sfr_z = take_log(sfr_z, volh, h0)
    mmol_z = take_log(mmol_z, volh, h0)
    mdust_z = take_log(mdust_z, volh, h0)
    sfr_tot = take_log(sfr_tot, volh, h0)
    mmol_tot = take_log(mmol_tot, volh, h0)
    mdust_tot = take_log(mdust_tot, volh, h0)

    mstar_mf[0,:,0] = xmf
    mstar_mf[1,:,0] = xmf
    mstar_mf[2,:,0] = xmf

    fname = 'Lagos19_quenched_smf.txt'
    np.savetxt(fname, mstar_mf[0,:], fmt='%.8e', delimiter=' ', newline='\n', header='', footer='', comments='# ')
    fname = 'Lagos19_starforming_smf.txt'
    np.savetxt(fname, mstar_mf[1,:], fmt='%.8e', delimiter=' ', newline='\n', header='', footer='', comments='# ')
    fname = 'Lagos19_allgals_smf.txt'
    np.savetxt(fname, mstar_mf[2,:], fmt='%.8e', delimiter=' ', newline='\n', header='', footer='', comments='# ')

    plot_sfr_contribution(plt, outdir, obsdir, sfr_z, sfr_tot, mmol_z, mmol_tot, mdust_z, mdust_tot, h0)
    plot_occupation_massive_halos(plt, outdir, obsdir, frac_mvir_occupation)
    plot_smf(plt, outdir, obsdir, mstar_mf)

if __name__ == '__main__':
    main(*common.parse_args())
