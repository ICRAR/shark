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

import numpy as np

import common
import utilities_statistics as us

mlow = 10
mupp = 15
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0


dm2 = 0.4
mbins2 = np.arange(mlow,mupp,dm2)
xmf2 = mbins2 + dm2/2.0



def prepare_data(hdf5_data, index, massgal, massbar, massbar_inside, massgal_morph, masshalo_massivegals, redshift, massgal_witherror):

    Omegab = 0.0491
    OmegaM = 0.3121

    (h0, _, mdisk, mbulge, mBH, mgas, mgas_bulge, mhot,
     mreheated, mhalo, typeg, idtree, sfrd, sfrb) = hdf5_data

    print("number of unique halos", len(np.unique(idtree)))
    #select most massive halos
    centrals = np.where(typeg == 0)
    print("number of halos with types==0", len(mhalo[centrals]))


    thresh = [0.5, 0.75]

    ind = np.where((typeg <= 0) & (mdisk+mbulge > 0))
    massgal[index,:] = us.wmedians(x=np.log10(mhalo[ind]) - np.log10(float(h0)),
                                   y=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                   xbins=xmf)

    ind = np.where((typeg <= 0) & (mdisk+mbulge > 0))
    massgal_witherror[index,:] = us.wmedians(x=np.log10(mhalo[ind]) - np.log10(float(h0)),
                                   y=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)) + np.random.normal(0.0, 0.3, len(mdisk[ind])),
                                   xbins=xmf)

    massbar[index,:] = us.wmedians(x=np.log10(mhalo[ind]) - np.log10(float(h0)),
                                   y=np.log10(mdisk[ind]+mbulge[ind]+mBH[ind]+mgas[ind]+mgas_bulge[ind]+mhot[ind]+mreheated[ind]) - np.log10(mhalo[ind]) - np.log10(Omegab/OmegaM),
                                   xbins=xmf)
    massbar_inside[index,:] = us.wmedians(x=np.log10(mhalo[ind]) - np.log10(float(h0)),
                                   y=np.log10(mdisk[ind]+mbulge[ind]+mBH[ind]+mgas[ind]+mgas_bulge[ind]+mhot[ind]) - np.log10(mhalo[ind]) - np.log10(Omegab/OmegaM),
                                   xbins=xmf)



    for i, c in enumerate(thresh):    
        ind = np.where((typeg <= 0) & (mdisk/(mdisk+mbulge) >= thresh[i]) & (mdisk+mbulge > 1e8))
        massgal_morph[i,0,index,:] = us.wmedians_2575(x=np.log10(mhalo[ind]) - np.log10(float(h0)),
                                       y=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                       xbins=xmf2, nmin=50)
    
        ind = np.where((typeg <= 0) & (mdisk/(mdisk+mbulge) < thresh[i]) & (mdisk+mbulge > 1e8))
        massgal_morph[i,1,index,:] = us.wmedians_2575(x=np.log10(mhalo[ind]) - np.log10(float(h0)),
                                       y=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                       xbins=xmf2,  nmin=50)
 
    ssfr_thresh = -1 + 0.5* redshift
    if (ssfr_thresh > 0): 
        ssfr_thresh = 0

    ran_err = np.random.normal(0.0, 0.3, len(mdisk))
    ran_err2 = np.random.normal(0.0, 0.3, len(mdisk))
    mnew = np.log10((mdisk+mbulge)/h0) + ran_err
    sfrnew = np.log10((sfrd+sfrb)/h0) + ran_err2

    ind = np.where(np.log10((mdisk+mbulge)/h0) > 10.7)
    H, _ = np.histogram(np.log10(mhalo[ind])- np.log10(float(h0)), bins=np.append(mbins,mupp))
    masshalo_massivegals[index,0,:] = masshalo_massivegals[index,0,:] + H
 
    ind = np.where(mnew > 10.7)
    H, _ = np.histogram(np.log10(mhalo[ind])- np.log10(float(h0)), bins=np.append(mbins,mupp))
    masshalo_massivegals[index,1,:] = masshalo_massivegals[index,1,:] + H

    ind = np.where((np.log10((mdisk+mbulge)/h0) > 10.7) & ((sfrd+sfrb)/(mdisk+mbulge) < 10**ssfr_thresh))
    H, _ = np.histogram(np.log10(mhalo[ind])- np.log10(float(h0)), bins=np.append(mbins,mupp))
    masshalo_massivegals[index,2,:] = masshalo_massivegals[index,2,:] + H
 
    ind = np.where((mnew > 10.7) & (sfrnew - mnew < ssfr_thresh))
    H, _ = np.histogram(np.log10(mhalo[ind])- np.log10(float(h0)), bins=np.append(mbins,mupp))
    masshalo_massivegals[index,3,:] = masshalo_massivegals[index,3,:] + H


    return thresh

def plot_SMHM_z(plt, outdir, zlist, massgal, obsdir, massgal_morph, thresh, masshalo_massivegals, massgal_witherror):

    def plot_moster13(ax, z, labels, label):
        #Moster et al. (2013) abundance matching SMHM relation
        M10 = 11.590
        M11 = 1.195
        N10 = 0.0351
        N11 = -0.0247
        beta10 = 1.376
        beta11 = -0.826
        gamma10 = 0.608
        gamma11 = 0.329
        M1 = pow(10.0, M10 + M11 * z/(z+1))
        N = N10 + N11 * z/(z+1)
        beta = beta10 + beta11 * z/(z+1)
        gamma = gamma10 + gamma11 * z/(z+1)

        mh = pow(10.0,xmf)
        m = mh * 2*N * pow (( pow(mh/M1, -beta ) + pow(mh/M1, gamma)), -1)

        if not labels:
            ax.plot(xmf,np.log10(m),'r', linestyle='dashed', linewidth=3)
        else:
            ax.plot(xmf,np.log10(m),'r', linestyle='dashed', linewidth=3, label=label)

    def plot_berhoozi13(ax, z, labels, label):
        a = 1.0/(1.0+z)
        nu = np.exp(-4*a*a)
        log_epsilon = -1.777 + (-0.006*(a-1)) * nu
        M1= 11.514 + ( - 1.793 * (a-1) - 0.251 * z) * nu
        alpha = -1.412 + 0.731 * nu * (a-1)
        delta = 3.508 + (2.608*(a-1)-0.043 * z) * nu
        gamma = 0.316 + (1.319*(a-1)+0.279 * z) * nu
        Min = xmf-M1
        fx = -np.log10(pow(10,alpha*Min)+1.0)+ delta * pow(np.log10(1+np.exp(Min)),gamma) / (1+np.exp(pow(10,-Min)))
        f = -0.3+ delta * pow(np.log10(2.0),gamma) / (1+np.exp(1))

        m = log_epsilon + M1 + fx - f

        if not labels:
            ax.plot(xmf,m, 'b', linestyle='dashdot',linewidth=3)
        else:
            ax.plot(xmf,m, 'b', linestyle='dashdot',linewidth=3, label=label)


    def plot_observations_kravtsov18(ax):

        mh, sm = common.load_observation(obsdir, 'SMHM/SatKinsAndClusters_Kravtsov18.dat', [0,1])
        ax.errorbar(mh, sm, xerr=0.2, yerr=0.2, color='purple', marker='s',  ls='None', label='Kravtsov+18')

    def plot_observations_romeo20(ax):

        mh, smmh = common.load_observation(obsdir, 'SMHM/Romeo20_SMHM.dat', [0,1])
        sm = np.log10(10**smmh * 10**mh)
        ax.plot(mh, sm, color='orange', marker='*',  ls='None', label='Romeo+20', fillstyle='none')

    def plot_observations_romeo20_morph(ax):

        mh, smmh = common.load_observation(obsdir, 'SMHM/Romeo20_SMHM.dat', [0,1])
        sm = np.log10(10**smmh * 10**mh)
        ax.plot(mh, sm, color='DarkTurquoise', marker='*',  ls='None', label='Romeo+20 (LTGs)', fillstyle='none')

        mh, smmh = common.load_observation(obsdir, 'SMHM/Romeo20_SMHM_ETGs.dat', [0,1])
        sm = np.log10(10**smmh * 10**mh)
        ax.plot(mh, sm, color='OrangeRed', marker='*',  ms=6, ls='None', label='Romeo+20 (ETGs)')

    def plot_observations_taylor20(ax):
        sm, sml, smh, hm, hml, hmh = common.load_observation(obsdir, 'SMHM/Taylor20.dat', [0,1,2,3,4,5])
        ax.errorbar(np.log10(hm*1e12), sm, xerr=[np.log10(hm)-np.log10(hml), np.log10(hmh)-np.log10(hm)], yerr=[sm-sml, smh-sm], color='Salmon', marker='d',  ls='None', label='Taylor+20')

    def plot_observations_kravtsov18_morph(ax):
        mh, sm, mhl, mhu = common.load_observation(obsdir, 'SMHM/LTGs_Kravtsov18.dat', [0,1,2,3])
        ind = np.where(mhl == 0)
        mhl[ind]=mh[ind] - 0.2
        mhu[ind]=mh[ind] + 0.2
        ax.errorbar(mh, sm, xerr=[mh-mhl, mhu-mh], yerr=0.2, color='DarkTurquoise', marker='s',  ls='None', label='Kravtsov+18 LTGs')
 
        mh, sm, mhl, mhu = common.load_observation(obsdir, 'SMHM/ETGs_Kravtsov18.dat', [0,1,2,3])
        ax.errorbar(mh, sm, xerr=[mh-mhl, mhu-mh], yerr=0.2, color='OrangeRed', marker='s',  ls='None', label='Kravtsov+18 ETGs')

    def plot_correa(ax):
        mh, sm, sml, smu = common.load_observation(obsdir, 'SMHM/LTGs_Correa19.dat', [0,1,2,3])
        ax.errorbar(mh, sm, yerr=[sm-sml, smu-sm], color='DarkTurquoise', marker='d',  ls='None', fillstyle='none', label='Correa+20 LTGs')

        mh, sm, sml, smu = common.load_observation(obsdir, 'SMHM/ETGs_Correa19.dat', [0,1,2,3])
        ax.errorbar(mh, sm, yerr=[sm-sml, smu-sm], color='OrangeRed', marker='d',  ls='None', fillstyle='none', label='Correa+20 ETGs')



    fig = plt.figure(figsize=(9.7,11.7))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    xmin, xmax, ymin, ymax = 10.5, 15, 7, 13
    xleg = xmin + 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)


    subplots = (321, 322, 323, 324, 325, 326)
    all_labels = (('Shark', 'Moster+13', 'Behroozi+13'), )

    for i, (z, subplot) in enumerate(zip(zlist, subplots)):
        labels = all_labels[0]
        
        # z=0 ##################################
        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

        ax.tick_params(labelsize=13)
        ax.text(xleg, yleg, 'z=%s' % str(z))

        #Predicted SMHM
        ind = np.where(massgal[i,0,:] != 0)
        xplot = xmf[ind]
        yplot = massgal[i,0,ind]
        errdn = massgal[i,1,ind]
        errup = massgal[i,2,ind]
   
        if not labels:
            ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
            ax.errorbar(xplot, yplot[0], color='k')
        else:
            ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
            ax.errorbar(xplot, yplot[0], color='k', label=labels[0])

        plot_moster13(ax, z, labels, labels[1])
        plot_berhoozi13(ax, z, labels, labels[2])

        if labels:
            common.prepare_legend(ax, ['r','b','k'], loc=4)


    common.savefig(outdir, fig, 'SMHM_z.pdf')

    fig = plt.figure(figsize=(9.7,11.7))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    xmin, xmax, ymin, ymax = 10.5, 15, 7, 13
    xleg = xmin + 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (321, 322, 323, 324, 325, 326)
    all_labels = (('Shark v2.0', 'v1.1 (L18)', 'Moster+13', 'Behroozi+13', 'EAGLE'), )


    def plot_l18(i, ax, label, lstyle='dotted', allgals = True):
        if(allgals):
           mh, sm = common.load_observation(obsdir, 'Models/SharkVariations/SMHM_Lagos18.dat', [0,i+1]) 
           ind = np.where(sm != 0)
           ax.plot(mh[ind], sm[ind], 'k', linestyle=lstyle, label=label)
        else:
           mh, smd, smb = common.load_observation(obsdir, 'Models/SharkVariations/SMHM_Lagos18_Morph.dat', [0,1,2])
           ind = np.where(smd > 0)
           ax.plot(mh[ind], smd[ind], 'b', linestyle=lstyle )
           ind = np.where(smb > 0)
           ax.plot(mh[ind], smb[ind], 'r', linestyle=lstyle)


    def plot_eagle(z, ax, label):
        mh, sm, red = common.load_observation(obsdir, 'Models/EAGLE/SMHM.dat', [0,1,4]) 
        ind = np.where(red == z)
        ax.plot(mh[ind], sm[ind], 'k', linestyle='dotted', label=label)

    z_eagle = [0, 0.5, 1.0, 2.0, 3.0, 4.0]

    for i, (z, subplot) in enumerate(zip(zlist, subplots)):
        labels = all_labels[0]
        
        # z=0 ##################################
        ax = fig.add_subplot(subplot)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

        ax.tick_params(labelsize=13)
        ax.text(xleg, yleg, 'z=%s' % str(z))

        #Predicted SMHM
        ind = np.where(massgal[i,0,:] != 0)
        xplot = xmf[ind]
        yplot = massgal[i,0,ind]
        errdn = massgal[i,1,ind]
        errup = massgal[i,2,ind]
  
        if not labels:
            ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='r', alpha=0.5, interpolate=True)
            ax.plot(xplot, yplot[0], color='darkred', linestyle='solid')
        else:
            ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='r', alpha=0.5,interpolate=True)
            ax.plot(xplot, yplot[0], color='darkred', linestyle='solid', label=labels[0])

        #ind = np.where(massgal_witherror[i,0,:] != 0)
        #xplot = xmf[ind]
        #yplot = massgal_witherror[i,0,ind]
        #if not labels:
        #    ax.plot(xplot, yplot[0], color='r', linestyle='dashed')
        #else:
        #    ax.plot(xplot, yplot[0], color='k', linestyle='dashed', label=labels[0] + '+0.3dex err')

        plot_l18(i, ax, labels[1], lstyle='dashed')
        #plot_eagle(z_eagle[i], ax, labels[4])
        #plot_moster13(ax, z, labels, labels[2])
        #plot_berhoozi13(ax, z, labels, labels[3])

        #if(i == 0):
        #   plot_observations_kravtsov18(ax)
        #   plot_observations_taylor20(ax)
        if labels:
            common.prepare_legend(ax, ['darkred','k','k','r','b'], loc=4)


    common.savefig(outdir, fig, 'SMHM_z_compL18.pdf')

    fig = plt.figure(figsize=(7,7))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    xmin, xmax, ymin, ymax = 10.5, 15, 7, 13
    xleg = xmin + 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (221, 222, 223, 224)
    all_labels = (('Shark v2.0', 'v1.1 (L18)', 'Moster+13', 'Behroozi+13', 'EAGLE'), )
    indices = [1,2,3,5]
  
    for j in range(0,len(indices)):
        labels = all_labels[0]
        z = zlist[indices[j]]
        i = indices[j]
        # z=0 ##################################
        ax = fig.add_subplot(subplots[j])
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

        ax.tick_params(labelsize=13)
        ax.text(xleg, yleg, 'z=%s' % str(z), fontsize=12)

        #Predicted SMHM
        ind = np.where(massgal[i,0,:] != 0)
        xplot = xmf[ind]
        yplot = massgal[i,0,ind]
        errdn = massgal[i,1,ind]
        errup = massgal[i,2,ind]
  
        if not labels:
            ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='r', alpha=0.5, interpolate=True)
            ax.plot(xplot, yplot[0], color='darkred', linestyle='solid')
        else:
            ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='r', alpha=0.5,interpolate=True)
            ax.plot(xplot, yplot[0], color='darkred', linestyle='solid', label=labels[0])

        #ind = np.where(massgal_witherror[i,0,:] != 0)
        #xplot = xmf[ind]
        #yplot = massgal_witherror[i,0,ind]
        #if not labels:
        #    ax.plot(xplot, yplot[0], color='r', linestyle='dashed')
        #else:
        #    ax.plot(xplot, yplot[0], color='k', linestyle='dashed', label=labels[0] + '+0.3dex err')

        plot_l18(i, ax, labels[1], lstyle='dashed')
        #plot_eagle(z_eagle[i], ax, labels[4])
        #plot_moster13(ax, z, labels, labels[2])
        #plot_berhoozi13(ax, z, labels, labels[3])

        #if(i == 0):
        #   plot_observations_kravtsov18(ax)
        #   plot_observations_taylor20(ax)
        if j == 0:
            common.prepare_legend(ax, ['darkred','k','k','r','b'], loc=4)

    plt.tight_layout()
    common.savefig(outdir, fig, 'SMHM_z_compL18_reduced.pdf')


    fig = plt.figure(figsize=(6,7))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    xleg = xmin + 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    z=0
    i=0
    # z=0 ##################################
    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    ax.tick_params(labelsize=13)
    ax.text(xleg, yleg, 'z=%s' % str(z))

    #Predicted SMHM
    ind = np.where(massgal[i,0,:] != 0)
    xplot = xmf[ind]
    yplot = massgal[i,0,ind]
    errdn = massgal[i,1,ind]
    errup = massgal[i,2,ind]
  
    if not labels:
        ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
        ax.plot(xplot, yplot[0], color='k', linestyle='solid')
    else:
        ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
        ax.plot(xplot, yplot[0], color='k', linestyle='solid', label=labels[0])

    plot_l18(i, ax, labels[1])
    plot_moster13(ax, z, labels, labels[2])
    plot_berhoozi13(ax, z, labels, labels[3])

    if(i == 0):
       plot_observations_kravtsov18(ax)
       plot_observations_taylor20(ax)
       plot_observations_romeo20(ax)
    if labels:
        common.prepare_legend(ax, ['k','k','r','b','orange','purple','salmon'], loc=4)

    ax = fig.add_subplot(212)
    xmin, xmax, ymin, ymax = 11.5, 14.8, 9, 12
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    plot_correa(ax)
    plot_observations_kravtsov18_morph(ax)
    #plot_observations_romeo20_morph(ax)
    #Predicted SMHM
    #ind = np.where(massgal[i,0,:] != 0)
    #xplot = xmf[ind]
    #yplot = massgal[i,0,ind]
    #errdn = massgal[i,1,ind]
    #errup = massgal[i,2,ind]
  
    #ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='grey', interpolate=True)
    #ax.plot(xplot, yplot[0], color='k', linestyle='solid', label='all')

    linestyles=['solid', 'dashed', 'dotted']

    i == 0
    ind = np.where(massgal_morph[i,0,0,0,:] != 0)
    xplot = xmf2[ind]
    yplot = massgal_morph[i,0,0,0,ind]
    errdn = massgal_morph[i,0,0,1,ind]
    errup = massgal_morph[i,0,0,2,ind]
    print(yplot[0]) 
    ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='blue', alpha=0.15, interpolate=True)
    ax.plot(xplot, yplot[0], color='b', linestyle=linestyles[i], label='D/T>0.5')
    #if i == 0:
    #   print("#Stellar-halo mass relation for B/T<0.5")
    #   for a,b in zip(xplot, yplot[0]):
    #       print(a,b)
   
    ind = np.where(massgal_morph[i,1,0,0,:] != 0)
    xplot = xmf2[ind]
    yplot = massgal_morph[i,1,0,0,ind]
    errdn = massgal_morph[i,1,0,1,ind]
    errup = massgal_morph[i,1,0,2,ind]
    
    ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='red', alpha=0.15, interpolate=True)
    ax.plot(xplot, yplot[0], color='r', linestyle=linestyles[i], label='D/T<0.5')
    #if i == 0:
    #   print("#Stellar-halo mass relation for B/T>0.5")
    #   for a,b in zip(xplot, yplot[0]):
    #       print(a,b)

    plot_l18(i, ax, 'v1.1', lstyle='dotted', allgals = False)
    common.prepare_legend(ax, ['b','r','DarkTurquoise','OrangeRed', 'DarkTurquoise','OrangeRed'], loc=4)

    plt.tight_layout()
    common.savefig(outdir, fig, 'SMHM_z0_compL18.pdf')


##############################################################################################################
################## halo mass distributions ###################################################################

    fig = plt.figure(figsize=(5,6))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit = "$\\rm PDF$"
    xleg = xmin + 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)
    xmin, xmax, ymin, ymax = 11.5, 14.8, 0, 4

    labels = ['z=0','z=0.5','z=1','z=2','z=3','z=4']

    cols = ['Lightblue','CornflowerBlue','DeepSkyBlue','DarkTurquoise','MediumBlue','MidnightBlue']
    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    ax.tick_params(labelsize=13)

    for i,l in enumerate(labels):
        #Predicted SMHM
        ind = np.where(masshalo_massivegals[i,0,:] != 0)
        xplot = xmf[ind]
        yplot = masshalo_massivegals[i,0,ind] / np.sum(masshalo_massivegals[i,0,ind] * dm)
        ax.plot(xplot, yplot[0], color=cols[i], linestyle='solid', label=labels[i])

        ind = np.where(masshalo_massivegals[i,1,:] != 0)
        xplot = xmf[ind]
        yplot = masshalo_massivegals[i,1,ind] / np.sum(masshalo_massivegals[i,1,ind] * dm)
        ax.plot(xplot, yplot[0], color=cols[i], linestyle='dashed')


    common.prepare_legend(ax, ['k','k','k','k','k','k'], loc=4)

    ax = fig.add_subplot(212)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    cols=['DarkOrange','LightSalmon','LightCoral','Crimson','Red','DarkRed']

    for i,l in enumerate(labels):
        #Predicted SMHM
        ind = np.where(masshalo_massivegals[i,2,:] != 0)
        xplot = xmf[ind]
        yplot = masshalo_massivegals[i,2,ind] / np.sum(masshalo_massivegals[i,2,ind] * dm)
        ax.plot(xplot, yplot[0], color=cols[i], linestyle='solid', label=labels[i])

        ind = np.where(masshalo_massivegals[i,3,:] != 0)
        xplot = xmf[ind]
        yplot = masshalo_massivegals[i,3,ind] / np.sum(masshalo_massivegals[i,3,ind] * dm)
        ax.plot(xplot, yplot[0], color=cols[i], linestyle='dashed')

    common.prepare_legend(ax, ['k','k','k','k','k','k'], loc=2)


    plt.tight_layout()
    common.savefig(outdir, fig, 'SMHM_z_massivegals.pdf')



def plot_BMHM_z(plt, outdir, massbar, massbar_inside):

    fig = plt.figure(figsize=(5,6))
    xtit = ""
    ytit = ""
    xmin, xmax, ymin, ymax = 10, 15, -1, 1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymin + 0.15 * (ymax - ymin)

    ax = fig.add_subplot(311)
    plt.subplots_adjust(left=0.17)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'z=0')

    #Predicted SMHM
    ind = np.where(massbar[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = massbar[0,0,ind]
    errdn = massbar[0,1,ind]
    errup = massbar[0,2,ind]

    ax.errorbar(xplot,yplot[0],color='k', label="all baryons")
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='+',markersize=2)

    ind = np.where(massbar_inside[0,0,:] != 0)
    xplot = xmf[ind]
    yplot = massbar_inside[0,0,ind]
    errdn = massbar_inside[0,1,ind]
    errup = massbar_inside[0,2,ind]

    ax.plot(xplot,yplot[0],color='b', linestyle='dotted', label="inside halos")
   
    xline = [10.0, 15.0]
    yline = [0.0, 0.0]
    ax.plot(xline,yline,'r', linestyle='dashed')

    common.prepare_legend(ax, ['b','k'], loc=2)

    # z=0.5 ##################################
    #ax = fig.add_subplot(222)
    #common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    #ax.text(xleg,yleg, 'z=0.5')

    #Predicted SMHM
    #ind = np.where(massbar[1,0,:] != 0)
    #xplot = xmf[ind]
    #yplot = massbar[1,0,ind]
    #errdn = massbar[1,1,ind]
    #errup = massbar[1,2,ind]

    #ax.errorbar(xplot,yplot[0],color='k', label="Shark")
    #ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='+',markersize=2)

    #ax.plot(xmf,xmf,'r', linestyle='dashed')


    # z=1 ##################################
    ax = fig.add_subplot(312)
    xtit = ""
    ytit = "$\\rm log_{10} (\\rm M_{\\rm bar}(\\Omega_{\\rm M}/\\Omega_{\\rm b})/\\rm M_{\\rm halo})$"

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg, yleg, 'z=1')

    # Predicted SMHM
    ind = np.where(massbar[2,0,:] != 0)
    xplot = xmf[ind]
    yplot = massbar[2,0,ind]
    errdn = massbar[2,1,ind]
    errup = massbar[2,2,ind]

    ax.errorbar(xplot,yplot[0],color='k', label="Shark")
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='+',markersize=2)

    ind = np.where(massbar_inside[2,0,:] != 0)
    xplot = xmf[ind]
    yplot = massbar_inside[2,0,ind]
    errdn = massbar_inside[2,1,ind]
    errup = massbar_inside[2,2,ind]

    ax.plot(xplot,yplot[0],color='b', linestyle='dotted')

    ax.plot(xline,yline,'r', linestyle='dashed')


    # z=1 ##################################
    ax = fig.add_subplot(313)
    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit = ""

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))
    ax.text(xleg,yleg, 'z=2')

    #Predicted SMHM
    ind = np.where(massbar[3,0,:] != 0)
    xplot = xmf[ind]
    yplot = massbar[3,0,ind]
    errdn = massbar[3,1,ind]
    errup = massbar[3,2,ind]

    ax.errorbar(xplot,yplot[0],color='k', label="Shark")
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='+',markersize=2)

    ind = np.where(massbar_inside[3,0,:] != 0)
    xplot = xmf[ind]
    yplot = massbar_inside[3,0,ind]
    errdn = massbar_inside[3,1,ind]
    errup = massbar_inside[3,2,ind]

    ax.plot(xplot,yplot[0],color='b', linestyle='dotted')

    ax.plot(xline,yline,'r', linestyle='dashed')

    common.savefig(outdir, fig, 'BMHM_z.pdf')


def main(modeldir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'm_bh', 'mgas_disk',
                           'mgas_bulge', 'mhot', 'mreheated', 'mvir_hosthalo',
                           'type', 'id_halo_tree', 'sfr_disk', 'sfr_burst')}

    zlist = (0, 0.5, 1, 2, 3, 4)
    snapshots = redshift_table[zlist]
    massgal = np.zeros(shape = (len(zlist), 3, len(xmf)))
    massgal_witherror = np.zeros(shape = (len(zlist), 3, len(xmf)))
    massgal_morph = np.zeros(shape = (4, 2, len(zlist), 3, len(xmf2)))
    massbar = np.zeros(shape = (len(zlist), 3, len(xmf)))
    massbar_inside =  np.zeros(shape = (len(zlist), 3, len(xmf)))
    masshalo_massivegals =  np.zeros(shape = (len(zlist), 4, len(xmf)))

    for idx, snapshot in enumerate(snapshots):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        thresh = prepare_data(hdf5_data, idx, massgal, massbar, massbar_inside, massgal_morph, masshalo_massivegals, zlist[idx], massgal_witherror)

    plot_SMHM_z(plt, outdir, zlist, massgal, obsdir, massgal_morph, thresh, masshalo_massivegals, massgal_witherror)
    plot_BMHM_z(plt, outdir, massbar, massbar_inside)

if __name__ == '__main__':
    main(*common.parse_args(requires_observations=True))
