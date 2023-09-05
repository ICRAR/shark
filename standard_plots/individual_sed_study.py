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
"""HMF plots"""

import numpy as np

import common


##################################
mlow = 5
mupp = 14
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3
c_light = 299792458.0 #m/s

def plot_individual_seds_z2(plt, outdir, obsdir, h0, SEDs_dust_z2, SEDs_app_z2, total_sfh_z2, disk_sfh_z2, sb_sfh_z2, gal_props_z2, LBT_z2):

    #wavelength in angstroms.
    file = obsdir+'/Models/Shark_SED_bands.dat'
    lambda_bands = np.loadtxt(file,usecols=[0],unpack=True)
    freq_bands   = c_light / (lambda_bands * 1e-10) #in Hz

    xtit="$\\rm log_{10}(\lambda/Ang\, (rest-frame))$"
    ytit="$\\rm log_{10}(\\nu f/ erg\, s^{-1}\, cm^{-2})$"

    xmin, xmax, ymin, ymax = 3.0, 7.0, -1, 6.5
    xleg = xmin + 0.1 * (xmax-xmin)
    yleg = ymax - 0.07 * (ymax-ymin)

    fig = plt.figure(figsize=(12,8.5))

    #subplots = (331, 332, 333, 334, 335, 336, 337, 338, 339)
    subplots = (231, 232, 233, 234, 235, 236)
    indices  = (0, 1, 2, 3, 4, 5) #, 6, 7, 8)
    sfrs     = (1000.0, 250.0, 100.0, 10.0, 5.0, 1.0)
    labelsfrs=  ('$\\rm SFR/M_{\odot} yr^{-1}= 500$', '$\\rm SFR/M_{\odot} yr^{-1}= 250$', '$\\rm SFR/M_{\odot} yr^{-1}= 100$', '$\\rm SFR/M_{\odot} yr^{-1}= 10$','$\\rm SFR/M_{\odot} yr^{-1}= 5$','$\\rm SFR/M_{\odot} yr^{-1}= 1$')

    for subplot, idx in zip(subplots, indices):

        ax = fig.add_subplot(subplot)
        ytitle = ' ' 
        xtitle = ' '
        if(idx == 0 or idx == 3):
           ytitle = ytit
        if(idx >= 3):
           xtitle = xtit
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitle, ytitle, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labelsfrs[idx])
        if(idx < 3):
           ax.text(3,6.6, '13.6')
           ax.text(4.7,6.7, 'LBT/Gyr')
           ax.text(6.9,6.6, '10.6')
        if(idx == 2 or idx == 5):
           ax.text(7.1, 4.5, '$\\rm log_{10}(SFR/M_{\odot} yr^{-1})$', rotation='vertical')

        ind = np.where((gal_props_z2[:,3] > sfrs[idx] - 0.2*sfrs[idx]) & (gal_props_z2[:,3] < sfrs[idx] + 0.2*sfrs[idx]) & (gal_props_z2[:,1] > 1e9))
        age_selec= gal_props_z2[ind,0]
        if(len(age_selec[0]) > 0):
           age_selec= gal_props_z2[ind,0]

           selecgal = 0 #len(age_selec[0])/2
           if idx > 3:
              selecgal = 1
           #choose first galaxy
           SEDs_selec   = SEDs_dust_z2[ind,4,:]
           SEDs_selec_b = SEDs_dust_z2[ind,2,:]
           SEDs_selec_d = SEDs_dust_z2[ind,3,:]
         
           tot_sfh_selec  = total_sfh_z2[ind,:]
           disk_sfh_selec = disk_sfh_z2[ind,:]
           sb_sfh_selec   = sb_sfh_z2[ind,:]

           SED_in   = pow(10.0, (SEDs_selec[0,selecgal] + 48.6)/(-2.5)) 
           SED_in_d = pow(10.0, (SEDs_selec_d[0,selecgal] + 48.6)/(-2.5)) 
           SED_in_b = pow(10.0, (SEDs_selec_b[0,selecgal] + 48.6)/(-2.5)) 

           yplot  = np.log10(SED_in * freq_bands)
           yplotd = np.log10(SED_in_d * freq_bands)
           yplotb = np.log10(SED_in_b * freq_bands)

           xplot = np.log10(lambda_bands) 
           ax.plot(xplot,yplot,'kx')
           ax.plot(xplot,yplotd,'bx')
           ax.plot(xplot,yplotb,'rx')

           sfh_in  = tot_sfh_selec[0,selecgal]
           sfh_ind = disk_sfh_selec[0,selecgal]
           sfh_inb = sb_sfh_selec[0,selecgal]

           lowsfr = np.where(sfh_in < 0.01)
           sfh_in[lowsfr] = 0.001
           lowsfr = np.where(sfh_ind < 0.01)
           sfh_ind[lowsfr] = 0.001
           lowsfr = np.where(sfh_inb < 0.01)
           sfh_inb[lowsfr] = 0.001

           maxLBT = max(LBT_z2)
           minLBT = min(LBT_z2)
           m = 4.0/(minLBT - maxLBT)
           b = 3.0 - m*maxLBT
           ax.plot(m*LBT_z2 + b, np.log10(sfh_in), 'k')
           ax.plot(m*LBT_z2 + b, np.log10(sfh_inb), 'r', linestyle='dashed')
           ax.plot(m*LBT_z2 + b, np.log10(sfh_ind), 'b', linestyle='dotted')

    common.savefig(outdir, fig, "SEDs_SFRs_z2.pdf")

    #wavelength in angstroms.
    file = obsdir+'/Models/Shark_SED_bands.dat'
    lambda_bands = np.loadtxt(file,usecols=[0],unpack=True)
    freq_bands   = c_light / (lambda_bands * 1e-10) #in Hz

    xtit="$\\rm log_{10}(\lambda/Ang\, (rest-frame))$"
    ytit="$\\rm log_{10}(\\nu f/ erg\, s^{-1}\, cm^{-2})$"

    xmin, xmax, ymin, ymax = 3.0, 7.0, -1, 6.8
    xleg = xmin + 0.1 * (xmax-xmin)
    yleg = ymax - 0.07 * (ymax-ymin)

    fig = plt.figure(figsize=(5,8))

    #subplots = (331, 332, 333, 334, 335, 336, 337, 338, 339)
    subplots = (311, 312, 313)
    indices  = (0, 1, 2) #, 6, 7, 8)
    sfrs     = (500.0, 100, 10.0)
    labelsfrs=  ('$\\rm SFR/M_{\odot} yr^{-1}= 1000$', '$\\rm SFR/M_{\odot} yr^{-1}= 100$', '$\\rm SFR/M_{\odot} yr^{-1}= 10$')

    for subplot, idx in zip(subplots, indices):

        ax = fig.add_subplot(subplot)
        xtitle = ' ' 
        ytitle = ytit
        if(idx >= 2):
           xtitle = xtit
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitle, ytitle, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labelsfrs[idx])
        if(idx== 0):
           ax.text(3,6.9, '13.6', fontsize=12)
           ax.text(4.7,7.1, 'LBT/Gyr', fontsize=12)
           ax.text(6.9,6.9, '10.6', fontsize=12)
        ax.text(7.1, 5.3, '$\\rm log_{10}(SFR/M_{\odot} yr^{-1})$', rotation='vertical', fontsize=12)

        ind = np.where((gal_props_z2[:,3] > sfrs[idx] - 0.2*sfrs[idx]) & (gal_props_z2[:,3] < sfrs[idx] + 0.2*sfrs[idx]) & (gal_props_z2[:,1] > 1e9))
        if(len(gal_props_z2[ind,3]) > 0):
           age_selec= gal_props_z2[ind,0]
           selecgal = len(age_selec)/2
           #if idx > 0:
           #   selecgal = 2
           #choose first galaxy
           SEDs_selec   = SEDs_dust_z2[ind,4,:]
           SEDs_selec_b = SEDs_dust_z2[ind,2,:]
           SEDs_selec_d = SEDs_dust_z2[ind,3,:]
        
           tot_sfh_selec  = total_sfh_z2[ind,:]
           disk_sfh_selec = disk_sfh_z2[ind,:]
           sb_sfh_selec   = sb_sfh_z2[ind,:]
           SED_in   = pow(10.0, (SEDs_selec[0,selecgal] + 48.6)/(-2.5)) 
           SED_in_d = pow(10.0, (SEDs_selec_d[0,selecgal] + 48.6)/(-2.5)) 
           SED_in_b = pow(10.0, (SEDs_selec_b[0,selecgal] + 48.6)/(-2.5)) 

           yplot  = np.log10(SED_in * freq_bands)
           yplotd = np.log10(SED_in_d * freq_bands)
           yplotb = np.log10(SED_in_b * freq_bands)

           xplot = np.log10(lambda_bands) 
           ax.plot(xplot,yplot,'ks', alpha=0.5)
           ax.plot(xplot,yplotd,'bd', alpha=0.5)
           ax.plot(xplot,yplotb,'rD', alpha=0.5)

           sfh_in  = tot_sfh_selec[0,selecgal]
           sfh_ind = disk_sfh_selec[0,selecgal]
           sfh_inb = sb_sfh_selec[0,selecgal]

           lowsfr = np.where(sfh_in < 0.01)
           sfh_in[lowsfr] = 0.001
           lowsfr = np.where(sfh_ind < 0.01)
           sfh_ind[lowsfr] = 0.001
           lowsfr = np.where(sfh_inb < 0.01)
           sfh_inb[lowsfr] = 0.001

           maxLBT = max(LBT_z2)
           minLBT = min(LBT_z2)
           m = 4.0/(minLBT - maxLBT)
           b = 3.0 - m*maxLBT
           ax.plot(m*LBT_z2 + b, np.log10(sfh_in), 'k')
           ax.plot(m*LBT_z2 + b, np.log10(sfh_inb), 'r', linestyle='dashed')
           ax.plot(m*LBT_z2 + b, np.log10(sfh_ind), 'b', linestyle='dotted')

    common.savefig(outdir, fig, "SEDs_SFRs_z2_threepanel.pdf")

    #wavelength in angstroms.
    xmin, xmax, ymin, ymax = 3.0, 7.0, -1, 6.1
    xleg = xmin + 0.1 * (xmax-xmin)
    yleg = ymax - 0.07 * (ymax-ymin)

    fig = plt.figure(figsize=(6,4.5))
    sfrs     = (100.0, 100, 10.0)
    colors   = ('Indigo','purple','Navy','DarkTurquoise', 'Aquamarine', 'Green','PaleGreen','GreenYellow','Gold','Yellow','Orange','OrangeRed','red','DarkRed','FireBrick','Crimson','IndianRed','LightCoral','Maroon','brown','Sienna','SaddleBrown','Chocolate','Peru','DarkGoldenrod','Goldenrod','SandyBrown')
    linestyles = ('dotted', 'dashed', 'solid')
    greys = ('Black','DarkSlateGray','DimGray')
    idx = 0
    ax = fig.add_subplot(111)
    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    left, bottom, width, height = 0.23, 0.25, 0.4, 0.2
    ax2 = fig.add_axes([left, bottom, width, height])

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(2, 2, 1, 1))
    plt.subplots_adjust(bottom = 0.2)
    #ax.text(xleg,yleg, '$\\rm SFR>250 M_{\\odot} yr^{-1}$')
    #ax.text(3,6.15, '13.6', fontsize=12)
    #ax.text(4.7,6.3, 'LBT/Gyr', fontsize=12)
    #ax.text(6.9,6.15, '10.6', fontsize=12)
    #ax.text(7.1, 4.3, '$\\rm log_{10}(SFR/M_{\odot} yr^{-1})$', rotation='vertical', fontsize=12)

    ind = np.where((gal_props_z2[:,3] > 255) & (gal_props_z2[:,1] > 1e9))
    age_selec= gal_props_z2[ind,0]
    SEDs_selec   = SEDs_dust_z2[ind,4,:]
    SEDs_app_selec = SEDs_app_z2[ind,4,:]
    print SEDs_app_selec
    tot_sfh_selec  = total_sfh_z2[ind,:]
    ssfrs = np.log10(gal_props_z2[ind,3]/gal_props_z2[ind,1]*1e9)
    ngals = 3
    for g in range(0,ngals):
        SED_in   = pow(10.0, (SEDs_selec[0,g] + 48.6)/(-2.5))
        yplot  = np.log10(SED_in * freq_bands)
        xplot = np.log10(lambda_bands) 
        ax.plot(xplot,yplot,color=greys[g],linestyle=linestyles[g], linewidth=1, label='$\\rm log_{10}(sSFR/Gyr^{-1})$=%.2f' % ssfrs[0,g]) 
        for xi,yi,c in zip(xplot,yplot,colors):
            ax.plot(xi,yi, 'd', markersize=4, color=c, alpha=0.5)
        sfh_in  = tot_sfh_selec[0,g]
        lowsfr = np.where(sfh_in < 0.01)
        sfh_in[lowsfr] = 0.001
 
        #maxLBT = max(LBT_z2)
        #minLBT = min(LBT_z2)
        #m = 4.0/(minLBT - maxLBT)
        #b = 3.0 - m*maxLBT
        ax2.plot(LBT_z2, np.log10(sfh_in), color=greys[g], linestyle=linestyles[g])

    common.prepare_legend(ax, greys, loc='center left')
    common.savefig(outdir, fig, "SEDs_SFRs_z2_starbursts.pdf")


def plot_individual_seds(plt, outdir, obsdir, h0, SEDs_dust_z0, SEDs_nodust_z0, total_sfh_z0, gal_props_z0, LBT):

    #wavelength in angstroms.
    file = obsdir+'/Models/Shark_SED_bands.dat'
    lambda_bands = np.loadtxt(file,usecols=[0],unpack=True)
    freq_bands   = c_light / (lambda_bands * 1e-10) #in Hz

    xtit="$\\rm log_{10}(\lambda/Ang\, (rest-frame))$"
    ytit="$\\rm log_{10}(\\nu f/ erg\, s^{-1}\, cm^{-2})$"

    xmin, xmax, ymin, ymax = 3.0, 7.0, -1, 4.5
    xleg = xmin + 0.1 * (xmax-xmin)
    yleg = ymax - 0.07 * (ymax-ymin)

    fig = plt.figure(figsize=(12,12))

    subplots = (331, 332, 333, 334, 335, 336, 337, 338, 339)
    indices  = (0, 1, 2, 3, 4, 5, 6, 7, 8)
    ages= (12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 5.0, 4.0, 3.0)
    labelages = ('12Gyr', '11Gyr', '10Gyr', '9Gyr','8Gyr','7Gyr','5Gyr','4Gyr','3Gyr')
    colors = ('Indigo','purple','Navy','DarkTurquoise', 'Aquamarine', 'Green','PaleGreen','GreenYellow','Gold','Yellow','Orange','OrangeRed','red','DarkRed','FireBrick','Crimson','IndianRed','LightCoral','Maroon','brown','Sienna','SaddleBrown','Chocolate','Peru','DarkGoldenrod','Goldenrod','SandyBrown')

    for subplot, idx in zip(subplots, indices):

        ax = fig.add_subplot(subplot)
        # These are in unitless percentages of the figure size. (0,0 is bottom left)
        left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]
        ax2 = fig.add_axes([left, bottom, width, height])

        ytitle = ' ' 
        xtitle = ' '
        if(idx == 0 or idx == 3 or idx == 6):
           ytitle = ytit
        if(idx >= 6):
           xtitle = xtit
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitle, ytitle, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labelages[idx])
        #if(idx < 3):
        #   ax.text(3,4.6, '13.6')
        #   ax.text(4.7,4.7, 'LBT/Gyr')
        #   ax.text(6.9,4.6, '0')
        #if(idx == 2 or idx == 5 or idx == 8):
        #   ax.text(7.1,-1.05,'-1')
        #   ax.text(7.1,3.95,'4')
        #   ax.text(7.1,3, '$\\rm log_{10}(SFR/M_{\odot} yr^{-1})$', rotation='vertical')

        ind = np.where((gal_props_z0[:,0] > ages[idx]-0.3) & (gal_props_z0[:,0] < ages[idx]+0.3) & (gal_props_z0[:,1] > 1e9))
        if(len(gal_props_z0[ind,0]) > 0):
           age_selec= gal_props_z0[ind,0]
           selecgal = 0 #len(age_selec[0])/2
           if idx == 3:
              selecgal = 1
           #choose first galaxy
           SEDs_selec        = SEDs_dust_z0[ind,4,:]
           SEDs_nodust_selec = SEDs_nodust_z0[ind,4,:]
           tot_sfh_selec = total_sfh_z0[ind,:]
           age_selec     = gal_props_z0[ind,0]  
           SED_in = pow(10.0, (SEDs_selec[0,selecgal] + 48.6)/(-2.5)) 
           SED_nodust_in = pow(10.0, (SEDs_nodust_selec[0,selecgal] + 48.6)/(-2.5))

           yplot   = np.log10(SED_in * freq_bands)
           yplotnd = np.log10(SED_nodust_in * freq_bands)
           xplot = np.log10(lambda_bands) 
   
           for xi,yi,yind,c in zip(xplot,yplot,yplotnd,colors):
               ax.plot(xi,yi, 'd', markersize=6, color=c, alpha=0.5)
               ax.plot(xi,yind, 'd', markersize=3, color=c, alpha=0.5)
  
           sfh_in = tot_sfh_selec[0,selecgal]
           lowsfr = np.where(sfh_in < 0.01)
           sfh_in[lowsfr] = 0.001
           ax2.plot(LBT, np.log10(sfh_in), 'r')

    common.savefig(outdir, fig, "SEDs_ages.pdf")
    
    xtit="$\\rm LBT/Gyr$"
    ytit="$\\rm log_{10}(SFR/M_{\odot} yr^{-1})$"

    xmin, xmax, ymin, ymax = 0, 13.6, -2, 3
    xleg = xmin + 0.1 * (xmax-xmin)
    yleg = ymax - 0.07 * (ymax-ymin)

    fig = plt.figure(figsize=(12,12))
    ages= (11.8, 11.0, 10.0, 9.0, 8.0, 7.0, 5.0, 4.0, 3.0)
    colors = ('Indigo','purple','Navy','DarkTurquoise', 'Aquamarine', 'Green','PaleGreen','GreenYellow','Gold','Yellow','Orange','OrangeRed','red','DarkRed','FireBrick','Crimson','IndianRed','LightCoral','Maroon','brown','Sienna','SaddleBrown','Chocolate','Peru','DarkGoldenrod','Goldenrod','SandyBrown')

    for subplot, idx in zip(subplots, indices):

        ax = fig.add_subplot(subplot)
        ytitle = ' ' 
        xtitle = ' '
        if(idx == 0 or idx == 3):
           ytitle = ytit
        if(idx >= 3):
           xtitle = xtit
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitle, ytitle, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labelages[idx])

        ind = np.where((gal_props_z0[:,0] > ages[idx]-0.3) & (gal_props_z0[:,0] < ages[idx]+0.3) & (gal_props_z0[:,1] > 1e9))
        if(len(gal_props_z0[ind]) > 0):
           age_selec= gal_props_z0[ind,0]
           selecgal = len(age_selec[0])/2
           tot_sfh_selec = total_sfh_z0[ind,:]
           age_selec     = gal_props_z0[ind,0]
           typesg        = gal_props_z0[ind,4]
           numgals       = 5
           if(numgals > len(age_selec[0])):
              numgals = len(age_selec[0])
 
           for gal in range(0,numgals): 
               sfh_in = tot_sfh_selec[0,gal]
               lowsfr = np.where(sfh_in < 0.01)
               sfh_in[lowsfr] = 0.001
               if(typesg[0,gal] == 0):
                  ax.plot(LBT, np.log10(sfh_in), color=colors[gal], linewidth=1)
               else:
                  ax.plot(LBT, np.log10(sfh_in), color=colors[gal], linewidth=1, linestyle = 'dashed')


    common.savefig(outdir, fig, "SFHs_ages.pdf")
    
    xtit="$\\rm log_{10}(\lambda/Ang\, (rest-frame))$"
    ytit="$\\rm log_{10}(\\nu f/ erg\, s^{-1}\, cm^{-2})$"

    xmin, xmax, ymin, ymax = 3.0, 7.0, -1, 4.5
    xleg = xmin + 0.1 * (xmax-xmin)
    yleg = ymax - 0.07 * (ymax-ymin)

    fig = plt.figure(figsize=(5,8.5))

    subplots = (311, 312, 313)
    indices  = (0, 1, 2)
    ages= (11.0, 8.0, 3.0)
    labelages = ('11Gyr', '8Gyr','3Gyr')
    colors = ('Indigo','purple','Navy','MediumBlue','Green','MediumAquamarine','LightGreen','YellowGreen','Gold','Orange','Coral','OrangeRed','red','DarkRed','FireBrick','Crimson','IndianRed','LightCoral','Maroon','brown','Sienna','SaddleBrown','Chocolate','Peru','DarkGoldenrod','Goldenrod','SandyBrown')

    botooms = [0.685,0.415,0.145]
    for subplot, idx in zip(subplots, indices):

        ax = fig.add_subplot(subplot)
        # These are in unitless percentages of the figure size. (0,0 is bottom left)
        left, bottom, width, height = 0.23, botooms[idx], 0.25, 0.07
        ax2 = fig.add_axes([left, bottom, width, height])

        xtitle = ' '
        ytitle = ytit
        if(idx == 2):
           xtitle = xtit
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitle, ytitle, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labelages[idx])

        ind = np.where((gal_props_z0[:,0] > ages[idx]-0.3) & (gal_props_z0[:,0] < ages[idx]+0.3) & (gal_props_z0[:,1] > 1e9))
        if(len(gal_props_z0[ind,0]) > 0):
           age_selec= gal_props_z0[ind,0]
           if(idx == 0 or idx == 2):
              selecgal = len(age_selec[0])/5
           else:
              selecgal = len(age_selec[0])/4
           #if idx == 3:
           #   selecgal = 1
           #choose first galaxy
           SEDs_selec        = SEDs_dust_z0[ind,4,:]
           SEDs_nodust_selec = SEDs_nodust_z0[ind,4,:]
           tot_sfh_selec = total_sfh_z0[ind,:]
           age_selec     = gal_props_z0[ind,0]  
           SED_in = pow(10.0, (SEDs_selec[0,selecgal] + 48.6)/(-2.5)) 
           SED_nodust_in = pow(10.0, (SEDs_nodust_selec[0,selecgal] + 48.6)/(-2.5))

           yplot   = np.log10(SED_in * freq_bands)
           yplotnd = np.log10(SED_nodust_in * freq_bands)
           xplot = np.log10(lambda_bands) 
   
           for xi,yi,yind,c in zip(xplot,yplot,yplotnd,colors):
               ax.plot(xi,yi, 'd', markersize=6, color=c, alpha=0.2)
               ax.plot(xi,yind, 'd', markersize=3, color=c)
  
           sfh_in = tot_sfh_selec[0,selecgal]
           lowsfr = np.where(sfh_in < 0.01)
           sfh_in[lowsfr] = 0.001
           ax2.plot(LBT, np.log10(sfh_in), 'r')

    common.savefig(outdir, fig, "SEDs_ages_threepanel.pdf")
    
    xtit="$\\rm LBT/Gyr$"
    ytit="$\\rm log_{10}(SFR/M_{\odot} yr^{-1})$"

    xmin, xmax, ymin, ymax = 0, 13.6, -2, 3
    xleg = xmin + 0.1 * (xmax-xmin)
    yleg = ymax - 0.07 * (ymax-ymin)

    fig = plt.figure(figsize=(5,8))
    ages= (11.0, 8.0, 3.0)
    colors = ('Indigo','purple','Navy','MediumBlue','Green','MediumAquamarine','LightGreen','YellowGreen','Gold','Orange','Coral','OrangeRed','red','DarkRed','FireBrick','Crimson','IndianRed','LightCoral','Maroon','brown','Sienna','SaddleBrown','Chocolate','Peru','DarkGoldenrod','Goldenrod','SandyBrown')

    for subplot, idx in zip(subplots, indices):

        ax = fig.add_subplot(subplot)
        ytitle = ytit
        xtitle = ' '
        if(idx >= 2):
           xtitle = xtit
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtitle, ytitle, locators=(2, 2, 1, 1))
        ax.text(xleg,yleg, labelages[idx])

        ind = np.where((gal_props_z0[:,0] > ages[idx]-0.3) & (gal_props_z0[:,0] < ages[idx]+0.3) & (gal_props_z0[:,1] > 1e9) & (gal_props_z0[:,4] <= 1))
        if(len(gal_props_z0[ind]) > 0):
           age_selec     = gal_props_z0[ind,0]
           tot_sfh_selec = total_sfh_z0[ind,:]
           age_selec     = gal_props_z0[ind,0]
           typesg        = gal_props_z0[ind,4]
           start         = len(age_selec[0])/2
           numgals       = 10
           if(numgals+start > len(age_selec[0])):
              numgals = len(age_selec[0])
              start = 0

           for gal in range(0,numgals):
               j = gal + start 
               sfh_in = tot_sfh_selec[0,j]
               lowsfr = np.where(sfh_in < 0.01)
               sfh_in[lowsfr] = 0.001
               if(typesg[0,j] == 0):
                  ax.plot(LBT, np.log10(sfh_in), color=colors[gal], linewidth=1)
               else:
                  ax.plot(LBT, np.log10(sfh_in), color=colors[gal], linewidth=1, linestyle = 'dashed')


    common.savefig(outdir, fig, "SFHs_ages_threepanel.pdf")

def prepare_data(hdf5_data, sfh, phot_data, phot_data_nod, phot_data_app, index, nbands):
   
    #star_formation_histories and SharkSED have the same number of galaxies in the same order, and so we can safely assume that to be the case.
    #to select the same galaxies in galaxies.hdf5 we need to ask for all of those that have a stellar mass > 0, and then assume that they are in the same order.

    (h0, _, mdisk, mbulge, mhalo, mshalo, typeg, age, 
     sfr_disk, sfr_burst, id_gal) = hdf5_data
   
    (bulge_diskins_hist, bulge_mergers_hist, disk_hist) = sfh

    #components:
    #(len(my_data), 2, 2, 5, nbands)
    #0: disk instability bulge
    #1: galaxy merger bulge
    #2: total bulge
    #3: disk
    #4: total
    #ignore last band which is the top-hat UV of high-z LFs.
    ind = np.where(mdisk + mbulge > 0)
    SEDs_dust   = np.zeros(shape = (len(mdisk[ind]), 5, nbands - 1))
    SEDs_nodust = np.zeros(shape = (len(mdisk[ind]), 5, nbands - 1))
    SEDs_app    = np.zeros(shape = (len(mdisk[ind]), 5, nbands - 1))
    print ("number of galaxies with mstar>0 %d" % len(mdisk[ind]))
    # print ("number of input galaxies in photometry file %d", len(phot_data[5]))

    p = 0
    for c in range(0,5):
        indust = phot_data[p]
        innodust = phot_data_nod[p]
        appdust = phot_data_app[p]
        for i in range(0,nbands-1):
            SEDs_dust[:,c,i] = indust[i,:]
            SEDs_nodust[:,c,i] = innodust[i,:]
            SEDs_app[:,c,i] = appdust[i,:]
        p = p + 1


    ngals       = len(SEDs_dust[:,0,0])
    nsnap       = len(bulge_diskins_hist[0,:])
    total_sfh = np.zeros(shape = (ngals, nsnap))
    sb_sfh    = np.zeros(shape = (ngals, nsnap))
    disk_sfh  = np.zeros(shape = (ngals, nsnap))
    gal_props = np.zeros(shape = (ngals, 5))

    for s in range(0,nsnap):
        total_sfh[:,s] = bulge_diskins_hist[:,s] + bulge_mergers_hist[:,s] + disk_hist[:,s] #in Msun/yr
        sb_sfh[:,s]    = bulge_diskins_hist[:,s] + bulge_mergers_hist[:,s]
        disk_sfh[:,s]  = disk_hist[:,s]

    ind = np.where((mdisk+mbulge) > 0)
    gal_props[:,0] = 13.6-age[ind]
    gal_props[:,1] = mdisk[ind] + mbulge[ind]
    gal_props[:,2] = mbulge[ind] / (mdisk[ind] + mbulge[ind])
    gal_props[:,3] = (sfr_burst[ind] + sfr_disk[ind])/1e9/h0
    gal_props[:,4] = typeg[ind]

    return (SEDs_dust, SEDs_nodust, SEDs_app, total_sfh, sb_sfh, disk_sfh, gal_props)

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    Variable_Ext = True
    file_hdf5_sed = "Shark-SED-eagle-rr14-steep.hdf5"

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                           'mvir_subhalo', 'type', 'mean_stellar_age', 
                           'sfr_disk', 'sfr_burst', 'id_galaxy')}

    sfh_fields = {'bulges_diskins': ('star_formation_rate_histories'),
                  'bulges_mergers': ('star_formation_rate_histories'),
                  'disks': ('star_formation_rate_histories')}

    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}
    fields_sed_nod = {'SED/ab_nodust': ('bulge_d','bulge_m','bulge_t','disk','total')}
    fields_sed_app = {'SED/ap_dust': ('bulge_d','bulge_m','bulge_t','disk','total')}

    z = (0, 2) #0.5, 1, 1.5, 2, 3)
    snapshots = redshift_table[z]

    # Create histogram
    for index, snapshot in enumerate(snapshots):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        sfh, delta_t, LBT = common.read_sfh(model_dir, snapshot, sfh_fields, subvols)
        seds = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed)
        seds_nod = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_nod, subvols, file_hdf5_sed)
        seds_app = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed_app, subvols, file_hdf5_sed)

        nbands = len(seds[0])
 
        (seds_dust, seds_nodust, sedsapp, total_sfh, sb_sfh, disk_sfh, gal_props) = prepare_data(hdf5_data, sfh, seds, seds_nod, seds_app, index, nbands)

        h0, volh = hdf5_data[0], hdf5_data[1]
        if(index == 0):
            SEDs_nodust_z0 = seds_nodust
            SEDs_dust_z0 = seds_dust
            total_sfh_z0 = total_sfh
            gal_props_z0 = gal_props
            LBT_z0 = LBT
            plot_individual_seds(plt, outdir, obsdir, h0, SEDs_dust_z0, SEDs_nodust_z0, total_sfh_z0, gal_props_z0, LBT_z0)

        if(index == 1):
            SEDs_dust_z2 = seds_dust
            SEDs_app_z2 = sedsapp
            total_sfh_z2 = total_sfh
            disk_sfh_z2  = disk_sfh
            sb_sfh_z2    = sb_sfh
            gal_props_z2 = gal_props
            LBT_z2 = LBT
            plot_individual_seds_z2(plt, outdir, obsdir, h0, SEDs_dust_z2, SEDs_app_z2, total_sfh_z2, disk_sfh_z2, sb_sfh_z2, gal_props_z2, LBT_z2)
  
if __name__ == '__main__':
    main(*common.parse_args())
