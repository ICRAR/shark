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
import utilities_statistics as us
import common
import functools


##################################
mlow = 8
mupp = 15
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

# Constants
GyrToYr = 1e9
zsun = 0.0189
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3
c_light = 299792458.0 #m/s

def plot_fraction_hydro(plt, outdir, obsdir, f_q, z):

    xtit="$\\rm log_{10}(M_{\\rm vir}/M_{\odot})$"
    ytit="$\\rm f_{\\rm hot-halo}$"

    xmin, xmax, ymin, ymax = 10, 15.0, -0.05, 1.05
    xleg = xmin + 0.1 * (xmax-xmin)
    yleg = ymax - 0.07 * (ymax-ymin)

    fig = plt.figure(figsize=(5,5))
    colors = ('Navy','Blue','RoyalBlue','SkyBlue','LightSalmon','IndianRed','Crimson','Red','DarkRed')

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(2, 2, 0.15, 0.15))
    #ax.text(xleg,yleg, labelages[idx])

    for i,j in enumerate(z):
        ind = np.where(f_q[i,:] !=0)
        xin = xmf[ind]
        yin = f_q[i,ind]
        ax.plot(xin, yin[0], color=colors[i], linewidth=3, label = 'z=%s' % str(j))

    common.prepare_legend(ax, ['k', 'k'], loc=2)
    common.savefig(outdir, fig, "hothalo_fraction_z.pdf")
    
def plot_individual_seds(plt, outdir, obsdir, h0, total_sfh_z0, gal_props_z0, LBT, redshift):

    xtit="$\\rm LBT/Gyr$"
    ytit="$\\rm log_{10}(SFR/M_{\odot} yr^{-1})$"

    xmin = 0
    if(redshift == '2'):
        xmin=10.35
    xmax, ymin, ymax = 13.6, -3.2, 3.5
    xleg = xmax + 0.025 * (xmax-xmin)
    yleg = ymax - 0.07 * (ymax-ymin)

    fig = plt.figure(figsize=(6.5,5))
    #mbins =  (11.0, 11.2, 11.4, 11.6, 11.8, 12.0, 12.2, 12.4, 12.6)
    mbins = (9,9.25,9.5,9.75,10,10.25,10.5,10.75,11,11.25,11.5, 11.75, 12.5)
    colors = ('Navy','Blue','RoyalBlue','SkyBlue','Teal','DarkTurquoise','Aquamarine','Yellow', 'Gold',  'Orange','OrangeRed', 'LightSalmon', 'Crimson', 'Red', 'DarkRed')

    ax = fig.add_subplot(111)
    if(redshift == '0'):
       common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(2, 2, 1, 1))
    if(redshift == '2'):
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.5, 0.5, 1, 1))
    ax.text(xleg,yleg,'$\\rm log_{10}(M_{\\star}/M_{\\odot})$', fontsize=12)
    ax.text(0.5,2.8,'Shark v2.0 (this work)', fontsize=13)
    #ax.text(0.5,2.8,'Shark v1.1 (L18)', fontsize=13)

    N_max = 25
    for j in range(0,len(mbins)-1):
        ind = np.where((gal_props_z0[:,1] > 10**mbins[j]) & (gal_props_z0[:,1] < 10**mbins[j+1]) & (gal_props_z0[:,4] == 0))
        if(len(gal_props_z0[ind]) > 0):
           age_selec= gal_props_z0[ind,0]
           tot_sfh_selec = total_sfh_z0[ind,:]
           tot_sfh_selec = tot_sfh_selec[0,:]
           age_selec     = gal_props_z0[ind,0]
           typesg        = gal_props_z0[ind,4]
           numgals = len(age_selec[0])
           if(numgals >= 10):
              print("Stellar mass bin", mbins[j], " has ", numgals, " galaxies") 
              SFH_med = np.zeros(shape = (3,len(LBT)))
              for snap in range(0,len(LBT)):
                  SFH_med[0,snap] = np.median(tot_sfh_selec[:,snap])
                  SFH_med[1,snap] = np.percentile(tot_sfh_selec[:,snap],25)
                  SFH_med[2,snap] = np.percentile(tot_sfh_selec[:,snap],75)
                  if(SFH_med[0,snap] < 0.001):
                     SFH_med[0,snap] = 0.001
                     SFH_med[1,snap] = 0.001 - 0.001*0.2
                     SFH_med[2,snap] = 0.001 + 0.001*0.2
              ax.fill_between(LBT,np.log10(SFH_med[1,:]),np.log10(SFH_med[2,:]), facecolor=colors[j], alpha=0.5, interpolate=True)
              ax.plot(LBT,np.log10(SFH_med[0,:]), color=colors[j], linewidth=3, label='%s' % str(mbins[j]+0.125))
              print("#SFH for stellar mass", str(mbins[j]+0.125))
              for a,b in zip(LBT,np.log10(SFH_med[0,:])):
                  print(a,b)
              #if(j == len(mbins)-4):
              #     for g in range(0,numgals):
              #         ax.plot(LBT, np.log10(tot_sfh_selec[g,:]), linewidth=1, color=colors[j])

           #nloop = numgals
           #if(numgals > N_max):
           #    nloop = N_max
           #for gal in range(0,nloop): 
           #    sfh_in = tot_sfh_selec[0,gal,:]
           #    print(sfh_in.shape)
           #    lowsfr = np.where(sfh_in < 0.01)
           #    sfh_in[lowsfr] = 0.001
           #    ax.plot(LBT, np.log10(sfh_in), color=colors[j], linewidth=1)


    common.prepare_legend(ax, colors, bbox_to_anchor=(0.98, 0.1))
    plt.tight_layout()
    common.savefig(outdir, fig, "SFHs_massivegalaxies_z"+redshift+".pdf")

def plot_age_stellar_mass(plt, outdir, obsdir, age_stellar_mass, age90_stellar_mass, age25_stellar_mass):


    fig = plt.figure(figsize=(5,4.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm t_{\\star,50}/Gyr$"
    xmin, xmax, ymin, ymax = 9, 12, 1, 14

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 2))
    mL18, ageL18, age84L18, age16L18 = np.loadtxt(obsdir + '/Models/SharkVariations/AgeMstar_Lagos18.dat', unpack = True, usecols = [0,1,2,3])

    ax.plot(mL18, ageL18,color='black',lw=2,label="Shark v1.1 (sm-w)")
    ax.fill_between(mL18, age84L18, age16L18, facecolor='k', alpha=0.5, interpolate=True)


    ind = np.where(age_stellar_mass[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = age_stellar_mass[0,0,ind]
        errdn = age_stellar_mass[0,1,ind]
        errup = age_stellar_mass[0,2,ind]
        ax.plot(xplot,yplot[0],color='red',lw=3.5,label="Shark v2.0 (sm-w)")
        ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='r', alpha=0.5, interpolate=True)

    mv19, agev19, agev19dn, agev19up =np.loadtxt(obsdir + '/stellarpop/age_mstar_vandensande19.dat', unpack = True, usecols = [0,1,2,3])
    ax.errorbar(mv19, 10**agev19, yerr=[10**agev19 - 10**agev19dn, 10**agev19up - 10**agev19], marker = 's', color='teal', label = 'van de Sande+19 (rband-w)')
    
    common.prepare_legend(ax, ['k', 'r', 'teal'], loc = 2)
    plt.tight_layout()
    common.savefig(outdir, fig, 'age_stellar_mass_z0.pdf')

    fig = plt.figure(figsize=(5,4.5))
    ytit = "$\\rm t_{\\star,90}/Gyr$"
    xmin, xmax, ymin, ymax = 9, 12, 0, 14

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 2))
    mL18, ageL18, age84L18, age16L18 = np.loadtxt(obsdir + '/Models/SharkVariations/AgeMstar_Lagos18.dat', unpack = True, usecols = [0,1,2,3])

    ax.plot(mL18, ageL18,color='black',lw=2,label="Shark v1.1 (sm-w)")
    ax.fill_between(mL18, age84L18, age16L18, facecolor='k', alpha=0.5, interpolate=True)


    ind = np.where(age90_stellar_mass[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = age90_stellar_mass[0,0,ind]
        errdn = age90_stellar_mass[0,1,ind]
        errup = age90_stellar_mass[0,2,ind]
        ax.plot(xplot,yplot[0],color='red',lw=3.5,label="Shark v2.0 (sm-w)")
        ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='r', alpha=0.5, interpolate=True)

    for a,b,c,d in zip(xplot,yplot[0],yplot[0]+errup[0],yplot[0]-errdn[0]):
        print(a,b,c,d)

    common.prepare_legend(ax, ['k', 'r'], loc = 2)
    plt.tight_layout()
    common.savefig(outdir, fig, 'age90_stellar_mass_z0.pdf')

    fig = plt.figure(figsize=(5,4.5))
    ytit = "$\\rm t_{\\star,25}/Gyr$"
    xmin, xmax, ymin, ymax = 9, 12, 0, 14

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 2))
    #mL18, ageL18, age84L18, age16L18 = np.loadtxt(obsdir + '/Models/SharkVariations/AgeMstar_Lagos18.dat', unpack = True, usecols = [0,1,2,3])

    #ax.plot(mL18, ageL18,color='black',lw=2,label="Shark v1.1 (sm-w)")
    #ax.fill_between(mL18, age84L18, age16L18, facecolor='k', alpha=0.5, interpolate=True)


    ind = np.where(age25_stellar_mass[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf[ind]
        yplot = age25_stellar_mass[0,0,ind]
        errdn = age25_stellar_mass[0,1,ind]
        errup = age25_stellar_mass[0,2,ind]
        ax.plot(xplot,yplot[0],color='red',lw=3.5,label="Shark v2.0 (sm-w)")
        ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='r', alpha=0.5, interpolate=True)

    for a,b,c,d in zip(xplot,yplot[0],yplot[0]+errup[0],yplot[0]-errdn[0]):
        print(a,b,c,d)

    common.prepare_legend(ax, ['k', 'r'], loc = 2)
    plt.tight_layout()
    common.savefig(outdir, fig, 'age25_stellar_mass_z0.pdf')



def prepare_data(hdf5_data, sfh, index, f_q, read_hydroeq, LBT, delta_t, age_stellar_mass, age90_stellar_mass, age25_stellar_mass):
   
    #star_formation_histories and SharkSED have the same number of galaxies in the same order, and so we can safely assume that to be the case.
    #to select the same galaxies in galaxies.hdf5 we need to ask for all of those that have a stellar mass > 0, and then assume that they are in the same order.
    bin_it   = functools.partial(us.wmedians, xbins=xmf)

    if(read_hydroeq):
       (h0, _, mdisk, mbulge, mhalo, mshalo, typeg, age, 
        sfr_disk, sfr_burst, id_gal, mbh, mhot, mreheated,
        on_hydrostatic_eq) = hdf5_data
    else:
        (h0, _, mdisk, mbulge, mhalo, mshalo, typeg, age,
        sfr_disk, sfr_burst, id_gal, mbh, mhot, mreheated) = hdf5_data


    if(read_hydroeq):
       for i,m in enumerate(xmf):
           ind = np.where((typeg == 0) & (on_hydrostatic_eq >= 0) & (mhalo/h0 > 10**(m-dm/2.0)) & (mhalo/h0 < 10**(m+dm/2.0)))
           n_all = len(typeg[ind])
           if(n_all > 9):
              ind = np.where((typeg == 0) & (on_hydrostatic_eq == 1) & (mhalo/h0 > 10**(m-dm/2.0)) & (mhalo/h0 < 10**(m+dm/2.0)))
              n_hydro = len(typeg[ind])
              f_q[index,i] = (n_hydro + 0.0) / (n_all + 0.0)
           #print(m, n_hydro, n_all)
   
   
       ind =np.where((typeg == 0) & (mhalo > 1e14))
       #print("massive halos", (on_hydrostatic_eq[ind]))

    (bulge_diskins_hist, bulge_mergers_hist, disk_hist) = sfh

    sfr_tot = (sfr_disk + sfr_burst)/1e9/h0

    #print(max(mhot/mhalo))
    #ind = np.where(mdisk + mbulge > 10**11.0)
    #print((sfr_burst[ind] + sfr_disk[ind])/1e9/h0)
    #print(np.log10(mbh[ind]), (mhot[ind]+mreheated[ind])/mhalo[ind], mdisk[ind]/(mdisk[ind]+mbulge[ind]))
  
    if(read_hydroeq):
       ind = np.where((mdisk + mbulge > 10**11.0) & (sfr_tot < 10) & (typeg == 0))
       #print(np.median(np.log10(mbh[ind])), np.median((mhot[ind]+mreheated[ind])/mhalo[ind]), np.median(mbulge[ind]/(mdisk[ind]+mbulge[ind])), np.median(np.log10(mbh[ind]/mbulge[ind])), np.median(np.log10(mhalo[ind])), np.median(on_hydrostatic_eq[ind]))
       Nlow = len(sfr_tot[ind])
   
       ind = np.where((mdisk + mbulge > 10**11.0) & (sfr_tot > 10) & (typeg == 0))
       #print(np.median(np.log10(mbh[ind])), np.median((mhot[ind]+mreheated[ind])/mhalo[ind]), np.median(mbulge[ind]/(mdisk[ind]+mbulge[ind])), np.median(np.log10(mbh[ind]/mbulge[ind])), np.median(np.log10(mhalo[ind])), np.median(on_hydrostatic_eq[ind]))
   
       Nhigh = len(sfr_tot[ind])
   
       #print("Fraction low over total", (Nlow+0.0)/(Nlow+Nhigh))
    #components:
    #(len(my_data), 2, 2, 5, nbands)
    #0: disk instability bulge
    #1: galaxy merger bulge
    #2: total bulge
    #3: disk
    #4: total
    #ignore last band which is the top-hat UV of high-z LFs.
    ind = np.where(mdisk + mbulge > 0)
    ngals       = len(mdisk[ind])
    nsnap       = len(bulge_diskins_hist[0,:])
    total_sfh = np.zeros(shape = (ngals, nsnap))
    sb_sfh    = np.zeros(shape = (ngals, nsnap))
    disk_sfh  = np.zeros(shape = (ngals, nsnap))
    gal_props = np.zeros(shape = (ngals, 6))

    gal_props[:,0] = 13.6-age[ind]
    gal_props[:,1] = mdisk[ind] + mbulge[ind]
    gal_props[:,2] = mbulge[ind] / (mdisk[ind] + mbulge[ind])
    gal_props[:,3] = (sfr_burst[ind] + sfr_disk[ind])/1e9/h0
    gal_props[:,4] = typeg[ind]
    
    for s in range(0,nsnap):
        total_sfh[:,s] = bulge_diskins_hist[:,s] + bulge_mergers_hist[:,s] + disk_hist[:,s] #in Msun/yr
        sb_sfh[:,s]    = bulge_diskins_hist[:,s] + bulge_mergers_hist[:,s]
        disk_sfh[:,s]  = disk_hist[:,s]

    tot_mass_formed = np.zeros(shape = (ngals))
    age_25 = np.zeros(shape = (ngals))
    age_50 = np.zeros(shape = (ngals))
    age_90 = np.zeros(shape = (ngals))

    for g in range(0,ngals):
        tot_mass_formed = sum(total_sfh[g,:] * delta_t * 1e9)
        mass_cum = 0
        s = 0
        while mass_cum < 0.9*tot_mass_formed:
            mass_cum += total_sfh[g,s] * delta_t[s] * 1e9
            if(mass_cum > 0.25*tot_mass_formed and age_25[g] == 0):
                age_25[g] = LBT[s]

            if(mass_cum > 0.5*tot_mass_formed and age_50[g] == 0):
                age_50[g] = LBT[s]
            s = s + 1
        age_90[g] = LBT[s-1]

    gal_props[:,5] = age_50

    age_stellar_mass[index,:] = bin_it(x=np.log10(gal_props[:,1]) - np.log10(float(h0)),
            y=age_50)
    age90_stellar_mass[index,:] = bin_it(x=np.log10(gal_props[:,1]) - np.log10(float(h0)),
            y=age_90)
    age25_stellar_mass[index,:] = bin_it(x=np.log10(gal_props[:,1]) - np.log10(float(h0)),
            y=age_25)

    return (total_sfh, sb_sfh, disk_sfh, gal_props)

def main(model_dir, outdir, redshift_table, subvols, obsdir):

    Variable_Ext = True
    file_hdf5_sed = "Shark-SED-eagle-rr14-steep.hdf5"
    read_hydroeq = False

    # Loop over redshift and subvolumes
    plt = common.load_matplotlib()
    if(read_hydroeq):
       fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                              'mvir_subhalo', 'type', 'mean_stellar_age', 
                              'sfr_disk', 'sfr_burst', 'id_galaxy', 'm_bh',
                              'mhot','mreheated','on_hydrostatic_eq')}
    else:
        fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mvir_hosthalo',
                              'mvir_subhalo', 'type', 'mean_stellar_age',
                              'sfr_disk', 'sfr_burst', 'id_galaxy', 'm_bh',
                              'mhot','mreheated')}


    sfh_fields = {'bulges_diskins': ('star_formation_rate_histories'),
                  'bulges_mergers': ('star_formation_rate_histories'),
                  'disks': ('star_formation_rate_histories')}

    z = (0, 2) #0.5, 1, 1.5, 2, 3)
    snapshots = redshift_table[z]

    f_q = np.zeros(shape = (len(z), len(xmf)))
    age_stellar_mass = np.zeros(shape = (len(z), 3, len(xmf)))
    age90_stellar_mass = np.zeros(shape = (len(z), 3, len(xmf)))
    age25_stellar_mass = np.zeros(shape = (len(z), 3, len(xmf)))

    # Create histogram
    for index, snapshot in enumerate(snapshots):

        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        sfh, delta_t, LBT = common.read_sfh(model_dir, snapshot, sfh_fields, subvols)

        (total_sfh, sb_sfh, disk_sfh, gal_props) = prepare_data(hdf5_data, sfh, index, f_q, read_hydroeq, LBT, delta_t, age_stellar_mass, age90_stellar_mass, age25_stellar_mass)

        h0, volh = hdf5_data[0], hdf5_data[1]
        total_sfh_z0 = total_sfh
        gal_props_z0 = gal_props
        LBT_z0 = LBT
        plot_individual_seds(plt, outdir, obsdir, h0, total_sfh_z0, gal_props_z0, LBT_z0, str(z[index]))

    plot_age_stellar_mass(plt, outdir, obsdir, age_stellar_mass, age90_stellar_mass, age25_stellar_mass)
    if(read_hydroeq):
       plot_fraction_hydro(plt, outdir, obsdir, f_q, z)
  
if __name__ == '__main__':
    main(*common.parse_args())
