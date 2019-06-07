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

import functools

import numpy as np
import os
import h5py

import common
import utilities_statistics as us

##################################
# Constants
mlow = 6.0
mupp = 12.0
dm = 0.2
mbins = np.arange(mlow, mupp, dm)
xmf = mbins + dm/2.0

zsun = 0.0189
MpctoKpc = 1e3

#spline fits to Trayford+2019 extinction paper

med_ar = [0.14414289, -1.01343542,  1.7168479]
perlow_ar = [0.12201443, -0.86945836,  1.47723534]
perhigh_ar = [ 0.14805331, -0.96745653,  1.48606172]

med_rv = [ 0.40421806, -2.38385783,  5.82676267]
perlow_rv = [  0.59000896,  -5.32212942,  14.93362796]
perhigh_rv = [  1.06070228,  -8.39519538,  21.01442891]


#read EAGLE tables
sdust_eaglet, taumed_eagle, taulow_eagle, tauhigh_eagle = common.load_observation('../data', 'Models/EAGLE/Tau5500-Trayford-EAGLE.dat', [0,1,2,3])
sdust_eaglem, mmed_eagle, mlow_eagle, mhigh_eagle = common.load_observation('../data/','Models/EAGLE/CFPowerLaw-Trayford-EAGLE.dat', [0,1,2,3])

print sdust_eaglet
m_med = np.zeros(shape = (2,len(sdust_eaglet)-1))
m_low = np.zeros(shape = (2,len(sdust_eaglet)-1))
m_hig = np.zeros(shape = (2,len(sdust_eaglet)-1))

for i in range(0,len(sdust_eaglet)-1):
    delta_dust = sdust_eaglet[i+1] - sdust_eaglet[i]
    m_med[0,i] = (taumed_eagle[i+1] - taumed_eagle[i] ) / delta_dust
    m_med[1,i] =  taumed_eagle[i+1]- m_med[0,i] * sdust_eaglet[i+1]
    print taumed_eagle[i+1], sdust_eaglet[i+1], m_med[0,i]*sdust_eaglet[i+1] + m_med[1,i]
    m_low[0,i] = (taulow_eagle[i+1] - taulow_eagle[i]) / delta_dust
    m_low[1,i] =  taulow_eagle[i+1]- m_low[0,i] * sdust_eaglet[i+1]
    m_hig[0,i] = (tauhigh_eagle[i+1] - tauhigh_eagle[i]) / delta_dust
    m_hig[1,i] =  tauhigh_eagle[i+1]- m_hig[0,i] * sdust_eaglet[i+1]

# define tau diffuse
def tau_diff (md, rd, hd, h0):

    tau    = np.zeros(shape = len(md))
    sigma  = np.zeros(shape = len(md))

    ind = np.where((md > 0) & (rd > 0) & (hd > 0))
    sigma[ind] = np.log10(md[ind]/h0 / (3.1416 * rd[ind]*MpctoKpc/h0 * hd[ind]*MpctoKpc/h0)) #in Msun/kpc^2
    # cap surface density of dust to physical values based on EAGLE
    sigma[ind] = np.clip(sigma[ind],min(sdust_eaglet), max(sdust_eaglet))

    #interpolate linearly in dust surface density going through the list of values in the EAGLE table
    for i in range(0,len(sdust_eaglet)-1):
        selecinrage = np.where((sigma >= sdust_eaglet[i]) & (sigma < sdust_eaglet[i+1]))
        tau[selecinrage] = m_med[0,i] * sigma[selecinrage] + m_med[1,i]
        var_gauss = ((m_low[0,i] * sigma[selecinrage] + m_low[1,i]) + (m_hig[0,i] * sigma[selecinrage] + m_hig[1,i]) ) * 0.5
        pert = np.random.randn(len(var_gauss)) * np.sqrt(var_gauss)
        tau[selecinrage] = tau[selecinrage] + pert

    # cap it to maximum and minimum values in EAGLE
    tau = np.clip(tau, 0.0, 3.5) 

    return (tau, sigma) 

# define clump tau
def tau_clump(mz,mg,tdiff):
    tau = np.zeros(shape = len(mz))
    ind = np.where((mz > 0) & (mg > 0))
    zgas = mz[ind]/mg[ind]
    zgas = np.clip(zgas, 1e-7, 2.5*zsun)
    tau[ind] = zgas/zsun

    # cap it to maximum and minimum values in EAGLE but also forcing the clump tau to be at least as high as the diffuse tau
    tau = np.clip(tau, tdiff, 3.5)
    return tau

def prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit):
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    xleg = xmax - 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)
    #ax.text(xleg, yleg, 'z=0')

def prepare_data(hdf5_data, index, tdiff, tcloud, sigmad_diff, sfr_rat, met_evo,  model_dir, snapshot, subvol, writeon):

    bin_it = functools.partial(us.wmedians, xbins=xmf)

    # Unpack data
    (h0, _, typeg, rgasd, rgasb, mHId, mH2d, mgasd, mHIb, mH2b, mgasb, mzd, mzb, mdisk, mbulge, sfrd, sfrb, idgal) = hdf5_data
    XH = 0.72
    h0log = np.log10(float(h0))

    # compute dust masses by assuming that a constant fraction of the metals are in dust.
    mdustd = 0.334 * mzd/h0
    mdustb = 0.334 * mzb/h0

    #random numbers from 0 to 1
    sin_inclination = np.random.rand(len(mdustd)) 
    inclination = np.arcsin(sin_inclination) * 180.0/3.1416 #degrees
    bd  = sin_inclination*(rgasd - rgasd/7.3)  +  rgasd/7.3 #scaleheight at r50


    (tau_dust_bulge, sigmab) = tau_diff(mdustb, 2.0*rgasb, 2.0*rgasb, h0)
    (tau_dust_disk, sigmad) = tau_diff(mdustd, 2.0*rgasd, 2.0*bd, h0)

    tau_clump_bulge = tau_clump(mzb, mgasb, tau_dust_bulge)
    tau_clump_disk  = tau_clump(mzd, mgasd, tau_dust_disk)

    mass = np.log10((mdisk +  mbulge)/h0)
    ind = np.where((mass >= 6) & (sfrd > 0))
    tdiff[index,0,:]  = bin_it(x=mass[ind], y=tau_dust_disk[ind])
    tcloud[index,0,:] = bin_it(x=mass[ind], y=tau_clump_disk[ind])
    sigmad_diff[index,0,:]  = bin_it(x=mass[ind], y=sigmad[ind])

    ind = np.where((mass >= 6) & (sfrb > 0))
    tdiff[index,1,:]  = bin_it(x=mass[ind], y=tau_dust_bulge[ind])
    tcloud[index,1,:] = bin_it(x=mass[ind], y=tau_clump_bulge[ind])
    sigmad_diff[index,1,:]  = bin_it(x=mass[ind], y=sigmab[ind])

    ind = np.where((mass >= 6) & (sfrd + sfrb > 0))
    sfr_rat[index,:] = bin_it(x=mass[ind], y=sfrb[ind]/(sfrd[ind]+sfrb[ind]))

    ind = np.where((mass >= 6) & (mgasd + mgasb > 0))
    met_evo[index,:] = bin_it(x=mass[ind], y=np.log10((mzd[ind]+mzb[ind])/(mgasd[ind]+mgasb[ind])/zsun))

    if(writeon):
        # will write the hdf5 files with the CO SLEDs and relevant quantities
        # will only write galaxies with mstar>0 as those are the ones being written in SFH.hdf5
        ind = np.where( (mdisk +  mbulge) > 0)
        file_to_write = os.path.join(model_dir, str(snapshot), str(subvol), 'extinction.hdf5')
        print ('Will write extinction to %s' % file_to_write)
        hf = h5py.File(file_to_write, 'w')
       
        hf.create_dataset('galaxies/tau_diff_disk', data=tau_dust_disk[ind])
        hf.create_dataset('galaxies/tau_diff_bulge', data=tau_dust_bulge[ind])
        hf.create_dataset('galaxies/tau_clump_disk', data=tau_clump_disk[ind])
        hf.create_dataset('galaxies/tau_clump_bulge', data=tau_clump_bulge[ind])
        hf.create_dataset('galaxies/id_galaxy', data=idgal[ind])
        hf.create_dataset('galaxies/inclination', data=inclination[ind])
        hf.close()

def plot_taus(plt, output_dir, tdiff, tcloud, sigmad_diff, sfr_rat, met_evo, zlist):

    #tau diffuse medium
    fig = plt.figure(figsize=(10,4.5))

    xmin, xmax, ymin, ymax = 6.0, 12.0, 0.0, 2.5
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\tau$"

    xleg = xmax - 0.6 * (xmax - xmin)
    yleg = ymax + 0.02 * (ymax - ymin)

    colors  = ('DarkRed','Salmon','Orange','YellowGreen','DarkTurquoise','MediumBlue','Purple')
    subplots = (131, 132, 133)
    labels = ('disk ISM','disk BC','bulge ISM')

    for s in range(0,len(subplots)):
        ax = fig.add_subplot(subplots[s])
        prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
        ax.text(xleg, yleg, labels[s])

        if(s == 0):
           j = 0
           for i in range(0,len(tcloud[:,j,0,0])):
               # Predicted relation
               ind = np.where(tdiff[i,j,0,:] > 0)
               xplot = xmf[ind]
               yplot = tdiff[i,j,0,ind]
               errdn = tdiff[i,j,1,ind]
               errup = tdiff[i,j,2,ind]
               ax.plot(xplot,yplot[0],color=colors[i])
               ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
               ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)
        if(s == 1):
           j = 0
           for i in range(0,len(tcloud[:,j,0,0])):
               # Predicted relation
               ind = np.where(tcloud[i,j,0,:] > 0)
               xplot = xmf[ind]
               yplot = tcloud[i,j,0,ind]
               errdn = tcloud[i,j,1,ind]
               errup = tcloud[i,j,2,ind]
               ax.plot(xplot,yplot[0],color=colors[i])
               ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
               ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)
        if(s == 2):
           j = 1
           for i in range(0,len(tcloud[:,j,0,0])):
               # Predicted relation
               ind = np.where(tdiff[i,j,0,:] > 0)
               xplot = xmf[ind]
               yplot = tdiff[i,j,0,ind]
               errdn = tdiff[i,j,1,ind]
               errup = tdiff[i,j,2,ind]
               ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
               ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
               ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)
        if(s == 0 or s == 2):
           x=[6,12]
           y=[0.5,0.5]
           ax.plot(x,y,linestyle='solid',color='k')
        else:
           x=[6,12]
           y=[1.5,1.5]
           ax.plot(x,y,linestyle='solid',color='k')

    common.prepare_legend(ax, colors, loc='upper left')
    common.savefig(output_dir, fig, "extinction_EAGLE_predictions.pdf")

    #tau clumps 
    fig = plt.figure(figsize=(7,4.5))

    xmin, xmax, ymin, ymax = 6.0, 12.0, 0.0, 2.5
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\tau_{\\rm clump}$"

    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (121, 122)
    labels = ('disk','bulge')

    for j in range(0,len(subplots)):
        ax = fig.add_subplot(subplots[j])
        prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
        ax.text(xleg, yleg, labels[j])

        for i in range(0,len(tcloud[:,j,0,0])):
            # Predicted relation
            ind = np.where(tcloud[i,j,0,:] > 0)
            xplot = xmf[ind]
            yplot = tcloud[i,j,0,ind]
            errdn = tcloud[i,j,1,ind]
            errup = tcloud[i,j,2,ind]
            ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
            ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
            ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)

        x=[6,12]
        y=[1.5,1.5]
        ax.plot(x,y,linestyle='solid',color='k')
    common.prepare_legend(ax, colors, loc='upper left')
    common.savefig(output_dir, fig, "extinction_clumps_EAGLE_predictions.pdf")


    #surface densities of dust
    fig = plt.figure(figsize=(7,4.5))

    xmin, xmax, ymin, ymax = 6.0, 12.0, 4, 8.1
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10} (\\Sigma_{\\rm dust}/M_{\odot} kpc^{-2})$"

    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (121, 122)
    labels = ('disk','bulge')

    for j in range(0,len(subplots)):
        ax = fig.add_subplot(subplots[j])
        prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
        ax.text(xleg, yleg, labels[j])

        for i in range(0,len(sigmad_diff[:,j,0,0])):
            # Predicted relation
            ind = np.where(sigmad_diff[i,j,0,:] > 0)
            xplot = xmf[ind]
            yplot = sigmad_diff[i,j,0,ind]
            errdn = sigmad_diff[i,j,1,ind]
            errup = sigmad_diff[i,j,2,ind]
            ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
            ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
            ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)

    common.prepare_legend(ax, colors, loc='upper left')
    common.savefig(output_dir, fig, "sigma_dust_predictions.pdf")

    #SFR ratio
    fig = plt.figure(figsize=(4.5,4.5))

    xmin, xmax, ymin, ymax = 6.0, 12.0, 0, 1
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm SFR_{\\rm burst}/SFR_{\\rm tot}$"

    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    for i in range(0,len(sfr_rat[:,0,0])):
        # Predicted relation
        ind = np.where(sfr_rat[i,0,:] > 0)
        xplot = xmf[ind]
        yplot = sfr_rat[i,0,ind]
        errdn = sfr_rat[i,1,ind]
        errup = sfr_rat[i,2,ind]
        ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)

    common.prepare_legend(ax, colors, loc='upper left')
    common.savefig(output_dir, fig, "sfr_ratio_predictions.pdf")

    #metallicity evolution
    fig = plt.figure(figsize=(4.5,4.5))

    xmin, xmax, ymin, ymax = 8.0, 12.0, -2, 0.8
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10} (Z_{\\rm gas}/Z_{\\odot})$"

    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(111)
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    plt.subplots_adjust(left=0.2)

    for i in range(0,len(sfr_rat[:,0,0])):
        # Predicted relation
        ind = np.where(met_evo[i,0,:] != 0)
        xplot = xmf[ind]
        yplot = met_evo[i,0,ind]
        errdn = met_evo[i,1,ind]
        errup = met_evo[i,2,ind]
        ax.plot(xplot,yplot[0],color=colors[i],label='z=%s' % str(zlist[i]))
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colors[i], alpha=0.5,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colors[i], alpha=0.5,interpolate=True)

    common.prepare_legend(ax, colors, loc='lower right')
    common.savefig(output_dir, fig, "metallicity_evo_predictions.pdf")

def main(model_dir, output_dir, redshift_table, subvols, obs_dir):

    zlist = (0, 0.5, 1, 2, 3, 4, 6)

    plt = common.load_matplotlib()
    fields = {'galaxies': ('type', 'rgas_disk', 'rgas_bulge', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge', 'mgas_metals_disk', 
                           'mgas_metals_bulge', 'mstars_disk', 'mstars_bulge','sfr_disk','sfr_burst','id_galaxy')}

    tau_diff = np.zeros(shape = (len(zlist), 2, 3, len(xmf)))
    tau_cloud = np.zeros(shape = (len(zlist), 2, 3, len(xmf)))
    sigmad_diff = np.zeros(shape = (len(zlist), 2, 3, len(xmf)))

    sfr_rat = np.zeros(shape = (len(zlist), 3, len(xmf)))
    met_evo = np.zeros(shape = (len(zlist), 3, len(xmf)))
    
    writeon = False

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        prepare_data(hdf5_data, index, tau_diff, tau_cloud, sigmad_diff, sfr_rat, met_evo, model_dir, snapshot, subvols, writeon)

    plot_taus(plt, output_dir, tau_diff, tau_cloud, sigmad_diff, sfr_rat, met_evo, zlist)

    writeon = True

    if(writeon):
       for index, snapshot in enumerate(redshift_table[zlist]):
           for subv in subvols:
               hdf5_data = common.read_data(model_dir, snapshot, fields, [subv])
               prepare_data(hdf5_data, index, tau_diff, tau_cloud, sigmad_diff, sfr_rat, met_evo, model_dir, snapshot, subv, writeon)

if __name__ == '__main__':
    main(*common.parse_args())
