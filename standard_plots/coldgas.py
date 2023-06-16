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

import common
import utilities_statistics as us

##################################
# Constants
mlow = 6.0
mupp = 12.0
dm = 0.5
mbins = np.arange(mlow, mupp, dm)
xmf = mbins + dm/2.0

slow = -1.5
supp = 2.5
ds = 0.2
sbins = np.arange(slow, supp, ds)
xsf = sbins + ds/2.0

sflow = -13.0
sfupp = -8.0
dsf = 0.5
sfbins = np.arange(sflow, sfupp, dsf)
xssfr = sfbins + dsf/2.0


def add_observations_to_plot(obsdir, fname, ax, marker, label, color='k', err_absolute=False):
    fname = '%s/Gas/%s' % (obsdir, fname)
    x, y, yerr_down, yerr_up = common.load_observation(obsdir, fname, (0, 1, 2, 3))
    common.errorbars(ax, x, y, yerr_down, yerr_up, color, marker, label=label, err_absolute=err_absolute)

def prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit):
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)
    xleg = xmax - 0.2 * (xmax-xmin)
    yleg = ymax - 0.1 * (ymax-ymin)
    #ax.text(xleg, yleg, 'z=0')

def prepare_data(index, redshift, hdf5_data):

    bin_it = functools.partial(us.wmedians, xbins=xmf)
    bin_it_ssfr = functools.partial(us.wmedians, xbins=xssfr)

    stack  = functools.partial(us.stacking, xbins=xmf)
    stack_sfr = functools.partial(us.stacking, xbins=xsf)
    stack_ssfr = functools.partial(us.stacking, xbins=xssfr)

    # Unpack data
    (h0, _, typeg, mdisk, mbulge, _, mbh, mHI, mH2, mgas,
     mHI_bulge, mH2_bulge, mgas_bulge, mvir, sfrd, 
     sfrb, id_halo_tree) = hdf5_data

    sfr = (sfrd + sfrb)/h0/1e9


    XH = 0.72
    h0log = np.log10(float(h0))

    n_typeg = len(typeg)
    morpho_type = np.zeros(shape = (n_typeg))
    morpho_type_stellar = np.zeros(shape = (n_typeg))
    ssfr = np.zeros(shape = (n_typeg))
    main_seq = np.zeros(shape = (n_typeg))
    type_g_conf = np.zeros(shape = (n_typeg))
    mhalo_conf = np.zeros(shape = (n_typeg))

    #we follow Campbell et al. (2015) who shows that the purity of satellites is constant with halo mass at a level of ~65%
    #while for centrals it goes from 90% at 1e12, 85% at 1e13, 80% at 1e14
    mhalo_bins=[0,1e12,1e13,1e14,1e16]
  
    CalculateHaloConfusion = False
    if(CalculateHaloConfusion == True):

       sat_pur = 0.65
       cen_put = [0.875, 0.875, 0.825, 0.7, 0.7]
       for j in range(0,len(mhalo_bins)-1):
           #select central galaxies 
           ind = np.where((typeg == 0)  & (mvir/h0 >= mhalo_bins[j]) & (mvir/h0 < mhalo_bins[j+1]))
           n_cen = len(typeg[ind])
           typein = np.zeros(shape = n_cen)
           #we now need to randomly select a fraction of centrals and reassigned them as satellites
           typein[:] = 0
           IDS = range(0,n_cen)
           selected = np.random.choice(IDS, size=int(np.floor((1 - cen_put[j+1])*n_cen)))
           typein[selected] = 1
           type_g_conf[ind] = typein
   
           #select satellites galaxies 
           ind = np.where((typeg > 0) & (mvir/h0 >= mhalo_bins[j]) & (mvir/h0 < mhalo_bins[j+1]))
           n_sat = len(typeg[ind])
           typein = np.zeros(shape = n_sat)
           #we now need to randomly select a fraction of satellites and reassigned them as centrals
           typein[:] = 1
           IDS = range(0,n_sat)
           selected = np.random.choice(IDS, size=int(np.floor((1-sat_pur)*n_sat)))
           typein[selected] = 0
           type_g_conf[ind] = typein
   
       #now we rank order the masses of halos and assign them to the new "centrals"
       #first take all central galaxies with a stellar mass >1e9Msun:
       ind = np.where((type_g_conf == 0) & ((mdisk + mbulge)/h0 > 1e9))
       n_most_mass_cens = len(type_g_conf[ind])
       ids_cen = id_halo_tree[ind]
       mstar_tot_group = np.zeros(shape = n_most_mass_cens)
       for i,g in enumerate(ids_cen):
           ind = np.where(id_halo_tree == g)
           mstar_tot_group[i] = sum(mdisk[ind] + mbulge[ind])
   
       #now select the n_most_mass_cens most massive halos
       ind = np.where(typeg == 0) 
       mvirin = mvir[ind]
       id_halo_tree_in = id_halo_tree[ind]
       mass_sorted_ids = np.argsort(1.0/mvirin)
       
       mhalo_most_massive = mvirin[mass_sorted_ids[0:n_most_mass_cens]]/h0
       ids_most_massive = id_halo_tree_in[mass_sorted_ids[0:n_most_mass_cens]]
   
       print(mhalo_most_massive, min(mhalo_most_massive), max(mhalo_most_massive))
   
       #sort total stellar mass of groups
       smass_sorted_ids = np.argsort(1.0/mstar_tot_group)
       ids_sm_sorted = ids_cen[smass_sorted_ids]
       for i,g in enumerate(ids_sm_sorted):
           ind = np.where(id_halo_tree == g)
           mhalo_conf[ind] = mhalo_most_massive[i]

    mhalo_conf = mvir/h0

    ssfr[:] = -15
    ind = np.where((mdisk + mbulge > 0) & (sfr > 0))
    ssfr[ind] = np.log10(sfr[ind]) - np.log10((mdisk[ind] + mbulge[ind])/h0)
    ind = np.where(ssfr + 9 > -1 + 0.5 * redshift)
    main_seq[ind] = 1
           
    ind = np.where((mbulge + mgas_bulge)/h0 > 1e12)
    print("median B/T massive gals", np.median(mbulge[ind]/(mbulge[ind] + mdisk[ind])))
    print("median SSFR massive gals", np.median(ssfr[ind]))
    print("median mBH massive gals", np.median(mbh[ind]))

    # Early-type galaxies criterion of Khochfar et al. (2011).
    ind = np.where((mbulge + mgas_bulge)/(mdisk + mbulge + mgas + mgas_bulge) > 0.5)
    morpho_type[ind] = 1.0

    # Early-type galaxies criterion based on stellar mass alone.
    ind = np.where((mbulge)/(mdisk + mbulge) > 0.5)
    morpho_type_stellar[ind] = 1.0

    mh2_gals = np.zeros(shape = (2, n_typeg))
    mh1_gals = np.zeros(shape = (2, n_typeg))
    mgas_gals = np.zeros(shape = (2, n_typeg))

    mh2_gals_ltg = np.zeros(shape = (2, n_typeg))
    mh1_gals_ltg = np.zeros(shape = (2, n_typeg))
    mgas_gals_ltg = np.zeros(shape = (2, n_typeg))

    mh2_gals_etg = np.zeros(shape = (2, n_typeg))
    mh1_gals_etg = np.zeros(shape = (2, n_typeg))
    mgas_gals_etg = np.zeros(shape = (2, n_typeg))

    mh1_relation_satellites_halos_stack = np.zeros(shape = (2, 5, 2, len(mbins))) 
    mh1_relation_ssfr_satellites_halos_stack = np.zeros(shape = (2, 5, 2, len(sfbins))) 
    mh1_relation_satellites_halos = np.zeros(shape = (2, 3, len(mbins)))
    mh1_relation_ssfr_satellites_halos = np.zeros(shape = (2, 3, len(sfbins)))

    #Gas scaling relation from stacking
    ind = np.where((mdisk + mbulge > 0) & (mHI_bulge + mHI > 0) & (main_seq == 1) & (sfr > 0))
    #print(mHI_bulge[ind] + mHI[ind])

    mh1_mstar_stack = stack(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=(mHI_bulge[ind] + mHI[ind])/h0)
    mh1_sfr_stack = stack_sfr(x=np.log10(sfr[ind]), y=(mHI_bulge[ind] + mHI[ind])/h0)
    fh1_mstar_stack = stack(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=(mHI_bulge[ind] + mHI[ind])/(mdisk[ind] + mbulge[ind]))
    fh1_sfr_stack = stack_sfr(x=np.log10(sfr[ind]), y=(mHI_bulge[ind] + mHI[ind])/(mdisk[ind] + mbulge[ind]))

    #print("HI stacking stellar mass at z=", redshift)    
    #for a,b,c in zip(xmf, mh1_mstar_stack, fh1_mstar_stack):
    #    print(a,b,c)
    #
    #print("HI stacking SFR at z=", redshift)  
    #for a,b,c in zip(xsf, mh1_sfr_stack, fh1_sfr_stack):
    #    print(a,b,c)



    #Gas scaling relations based on morphological criterion calculated using total baryon mass of disk and bulge
    ind = np.where((mdisk + mbulge > 0) & (mgas + mgas_bulge > 0) & (mH2 + mH2_bulge > 0) & (morpho_type == 0))
    # Data we'll use later
    mass = mdisk[ind] + mbulge[ind]
    mgas_gals_ltg[0,ind] = mh1_gals_ltg[0,ind] = mh2_gals_ltg[0,ind] = np.log10(mass) - h0log
    mgas_gals_ltg[1,ind] = np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / mass)
    mh1_gals_ltg[1,ind] = np.log10(XH * (mHI[ind] + mHI_bulge[ind]) / mass)
    mh2_gals_ltg[1,ind] = np.log10(XH * (mH2[ind] + mH2_bulge[ind]) / (mass))

    mgas_relation_ltg = bin_it(x=mgas_gals_ltg[0, ind], y=mgas_gals_ltg[1, ind])
    mh1_relation_ltg = bin_it(x=mh1_gals_ltg[0, ind], y=mh1_gals_ltg[1, ind])
    mh2_relation_ltg = bin_it(x=mh2_gals_ltg[0, ind], y=mh2_gals_ltg[1, ind])

    ind = np.where((mdisk + mbulge > 0) & (mgas + mgas_bulge > 0) & (mH2 + mH2_bulge > 0) & (morpho_type == 1))
    # Data we'll use later
    mass = mdisk[ind] + mbulge[ind]
    mgas_gals_etg[0,ind] = mh1_gals_etg[0,ind] = mh2_gals_etg[0,ind] = np.log10(mass) - h0log
    mgas_gals_etg[1,ind] = np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / mass)
    mh1_gals_etg[1,ind] = np.log10(XH * (mHI[ind] + mHI_bulge[ind]) / mass)
    mh2_gals_etg[1,ind] = np.log10(XH * (mH2[ind] + mH2_bulge[ind]) / (mass))

    mgas_relation_etg = bin_it(x=mgas_gals_etg[0, ind], y=mgas_gals_etg[1, ind])
    mh1_relation_etg = bin_it(x=mh1_gals_etg[0, ind], y=mh1_gals_etg[1, ind])
    mh2_relation_etg = bin_it(x=mh2_gals_etg[0, ind], y=mh2_gals_etg[1, ind])

    #Gas scaling relations based on morphological criterion calculated using stellar mass of disk and bulge
    ind = np.where((mdisk + mbulge > 0) & (mgas + mgas_bulge > 0) & (mH2 + mH2_bulge > 0) & (morpho_type_stellar == 0))
    # Data we'll use later
    mass = mdisk[ind] + mbulge[ind]
    mgas_gals_ltg[0,ind] = mh1_gals_ltg[0,ind] = mh2_gals_ltg[0,ind] = np.log10(mass) - h0log
    mgas_gals_ltg[1,ind] = np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / mass)
    mh1_gals_ltg[1,ind] = np.log10(XH * (mHI[ind] + mHI_bulge[ind]) / mass)
    mh2_gals_ltg[1,ind] = np.log10(XH * (mH2[ind] + mH2_bulge[ind]) / (mass))

    mgas_ms_relation_ltg = bin_it(x=mgas_gals_ltg[0, ind], y=mgas_gals_ltg[1, ind])
    mh1_ms_relation_ltg = bin_it(x=mh1_gals_ltg[0, ind], y=mh1_gals_ltg[1, ind])
    mh2_ms_relation_ltg = bin_it(x=mh2_gals_ltg[0, ind], y=mh2_gals_ltg[1, ind])

    ind = np.where((mdisk + mbulge > 0) & (mgas + mgas_bulge > 0) & (mH2 + mH2_bulge > 0) & (morpho_type_stellar == 1))
    # Data we'll use later
    mass = mdisk[ind] + mbulge[ind]
    mgas_gals_etg[0,ind] = mh1_gals_etg[0,ind] = mh2_gals_etg[0,ind] = np.log10(mass) - h0log
    mgas_gals_etg[1,ind] = np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / mass)
    mh1_gals_etg[1,ind] = np.log10(XH * (mHI[ind] + mHI_bulge[ind]) / mass)
    mh2_gals_etg[1,ind] = np.log10(XH * (mH2[ind] + mH2_bulge[ind]) / (mass))

    mgas_ms_relation_etg = bin_it(x=mgas_gals_etg[0, ind], y=mgas_gals_etg[1, ind])
    mh1_ms_relation_etg = bin_it(x=mh1_gals_etg[0, ind], y=mh1_gals_etg[1, ind])
    mh2_ms_relation_etg = bin_it(x=mh2_gals_etg[0, ind], y=mh2_gals_etg[1, ind])

    # Constrains
    ind = np.where((mdisk + mbulge > 0) & (mgas + mgas_bulge > 0) & (mH2 + mH2_bulge > 0))

    # Data we'll use later
    mass = mdisk[ind] + mbulge[ind]
    mgas_gals[0,ind] = mh1_gals[0,ind] = mh2_gals[0,ind] = np.log10(mass) - h0log
    mgas_gals[1,ind] = np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / mass)
    mh1_gals[1,ind] = np.log10(XH * (mHI[ind] + mHI_bulge[ind]) / mass)
    mh2_gals[1,ind] = np.log10(XH * (mH2[ind] + mH2_bulge[ind]) / (mass))

    # Binned relations
    mgas_relation = bin_it(x=mgas_gals[0, ind], y=mgas_gals[1, ind])
    mh1_relation = bin_it(x=mh1_gals[0, ind], y=mh1_gals[1, ind])
    mh2_relation = bin_it(x=mh2_gals[0, ind], y=mh2_gals[1, ind])
    mhr_relation = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - h0log,
                        y=np.log10((mH2[ind] + mH2_bulge[ind]) / (mHI[ind] + mHI_bulge[ind])))

    for i in range(0,len(mhalo_bins)-1):
        ind = np.where((mdisk + mbulge > 0) & (mvir/h0 >= mhalo_bins[i]) & (mvir/h0 < mhalo_bins[i+1]) & (typeg > 0) & ((mdisk + mbulge)/h0 > 1e9))
        mh1_relation_satellites_halos_stack[0,i+1,:] =  stack(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=(mHI[ind]+mHI_bulge[ind])/(mdisk[ind] + mbulge[ind]))
        mh1_relation_ssfr_satellites_halos_stack[0,i+1,:] =  stack_ssfr(x=ssfr[ind], y=(mHI[ind]+mHI_bulge[ind])/(mdisk[ind] + mbulge[ind]))
 
        ind = np.where((mdisk + mbulge > 0) & (mhalo_conf >= mhalo_bins[i]) & (mhalo_conf < mhalo_bins[i+1]) & (type_g_conf > 0) & ((mdisk + mbulge)/h0 > 1e9))
        mh1_relation_satellites_halos_stack[1,i+1,:] =  stack(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=(mHI[ind]+mHI_bulge[ind])/(mdisk[ind] + mbulge[ind]))
        mh1_relation_ssfr_satellites_halos_stack[1,i+1,:] =  stack_ssfr(x=ssfr[ind], y=(mHI[ind]+mHI_bulge[ind])/(mdisk[ind] + mbulge[ind]))
   
    ind = np.where((mdisk + mbulge > 0) & (typeg > 0) & ((mdisk + mbulge)/h0 > 1e9))
    mh1_relation_satellites_halos_stack[0,0,:] =  stack(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=(mHI[ind]+mHI_bulge[ind])/(mdisk[ind] + mbulge[ind]))
    mh1_relation_ssfr_satellites_halos_stack[0,0,:] =  stack_ssfr(x=ssfr[ind], y=(mHI[ind]+mHI_bulge[ind])/(mdisk[ind] + mbulge[ind]))
    mh1_relation_satellites_halos[0,:] = bin_it(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=(mHI[ind]+mHI_bulge[ind])/(mdisk[ind] + mbulge[ind]))
    mh1_relation_ssfr_satellites_halos[0,:] =  bin_it_ssfr(x=ssfr[ind], y=(mHI[ind]+mHI_bulge[ind])/(mdisk[ind] + mbulge[ind]))

    ind = np.where((mdisk + mbulge > 0) & (type_g_conf > 0) & ((mdisk + mbulge)/h0 > 1e9))
    mh1_relation_satellites_halos_stack[1,0,:] =  stack(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=(mHI[ind]+mHI_bulge[ind])/(mdisk[ind] + mbulge[ind]))
    mh1_relation_ssfr_satellites_halos_stack[1,0,:] =  stack_ssfr(x=ssfr[ind], y=(mHI[ind]+mHI_bulge[ind])/(mdisk[ind] + mbulge[ind]))
    mh1_relation_satellites_halos[1,:] = bin_it(x=np.log10((mdisk[ind] + mbulge[ind])/h0), y=(mHI[ind]+mHI_bulge[ind])/(mdisk[ind] + mbulge[ind]))
    mh1_relation_ssfr_satellites_halos[1,:] =  bin_it_ssfr(x=ssfr[ind], y=(mHI[ind]+mHI_bulge[ind])/(mdisk[ind] + mbulge[ind]))


    ind = np.where((mdisk+mbulge > 0) & (typeg == 0) & ((mdisk + mbulge)/h0 > 1e9))
    mass_central = mdisk[ind] + mbulge[ind]
    mgas_relation_cen = bin_it(x=np.log10(mass_central) - h0log,
                               y=np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / (mdisk[ind] + mbulge[ind] + XH*mgas[ind] + XH*mgas_bulge[ind])))
    mhr_relation_cen = bin_it(x=np.log10(mass_central) - h0log,
                              y=np.log10((mH2[ind] + mH2_bulge[ind]) / (mHI[ind] + mHI_bulge[ind])))

    ind = np.where((mdisk+mbulge > 0) & (typeg > 0))
    mass_sat = np.log10(mdisk[ind] + mbulge[ind]) - h0log
    mgas_relation_sat = bin_it(x=mass_sat,
                               y=np.log10(XH * (mgas[ind] + mgas_bulge[ind]) / (mdisk[ind] + mbulge[ind] + XH*mgas[ind] + XH*mgas_bulge[ind])))
    mhr_relation_sat = bin_it(x=mass_sat,
                              y=np.log10((mH2[ind] + mH2_bulge[ind]) / (mHI[ind] + mHI_bulge[ind])))

    return (mgas_relation, mgas_relation_cen, mgas_relation_sat,
            mh2_gals, mh1_gals, mgas_gals,
            mh2_relation, mh1_relation, mhr_relation, mhr_relation_cen, mhr_relation_sat, 
            mgas_relation_ltg, mh2_relation_ltg, mh1_relation_ltg, mgas_relation_etg, mh2_relation_etg, 
	    mh1_relation_etg, mgas_ms_relation_ltg, mh2_ms_relation_ltg, mh1_ms_relation_ltg, 
	    mgas_ms_relation_etg, mh2_ms_relation_etg, mh1_ms_relation_etg, mh1_relation_satellites_halos,
            mh1_relation_ssfr_satellites_halos, mh1_relation_satellites_halos_stack, mh1_relation_ssfr_satellites_halos_stack)

def plot_cold_gas_fraction(plt, output_dir, obs_dir, mgas_relation, mgas_relation_cen, mgas_relation_sat):

    ###################################
    #   Plots global mass densities
    fig = plt.figure(figsize=(5,4.5))

    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm cold}/M_{\\star})$"

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    prepare_ax(ax, 8, 12, -3, 0.1, xtit, ytit)

    #Predicted SMHM
    ind = np.where(mgas_relation[0,:] != 0)
    xplot = xmf[ind]
    yplot = mgas_relation[0,ind]
    errdn = mgas_relation[1,ind]
    errup = mgas_relation[2,ind]

    ax.errorbar(xplot,yplot[0],color='k', label="all galaxies")
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='+',markersize=2)

    ind = np.where(mgas_relation_cen[0,:] != 0)
    xplot = xmf[ind]
    yplot = mgas_relation_cen[0,ind]
    ax.errorbar(xplot,yplot[0],color='b',linestyle='dotted', label="centrals")

    ind = np.where(mgas_relation_sat[0,:] != 0)
    xplot = xmf[ind]
    yplot = mgas_relation_sat[0,ind]
    ax.errorbar(xplot,yplot[0],color='r',linestyle='dashed', label="satelites")

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    add_observations_to_plot(obs_dir, 'NeutralGasRatio_NonDetEQUpperLimits.dat', ax, 'v', "GASS+COLDGASS", color='grey')
    add_observations_to_plot(obs_dir, 'NeutralGasRatio_NonDetEQZero.dat', ax, '^', "GASS+COLDGASS", color='grey')

    common.prepare_legend(ax, ['k','b','r','grey','grey'])
    common.savefig(output_dir, fig, "cold_gas_fraction.pdf")

def plot_HI_stacking(plt, output_dir, obs_dir, mh1_relation_satellites_halos, mh1_relation_ssfr_satellites_halos):

    ###################################
    #   Plots global mass densities
    fig = plt.figure(figsize=(16,4.5))

    subp = [151, 152, 153, 154, 155]
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm HI}/M_{\\star})$"


    labels = ['all satellites', '$M_{\\rm halo}<10^{12}\\rm M_{\\odot}$', '$10^{12}\\rm M_{\\odot}<M_{\\rm halo}<10^{13}\\rm M_{\\odot}$', '$10^{13}\\rm M_{\\odot}<M_{\\rm halo}<10^{14}\\rm M_{\\odot}$', '$M_{\\rm halo}>10^{14}\\rm M_{\\odot}$']

    def plot_obs_brown17(ax, bin_mass=0, mass = True, lab=True, delta=False, color='Salmon'):
        if(mass == True):
           x, y, ydown, yup, binm = common.load_observation(obs_dir, 'Gas/HIMstar_Brown17.dat', (0, 1, 2, 3, 4))
        else:
           x, y, ydown, yup, binm = common.load_observation(obs_dir, 'Gas/HISSFR_Brown17.dat', (0, 1, 2, 3, 4))

        #define the relation for the whole sample
        ind = np.where(binm == 0)
        xin_all = x[ind]
        yin_all = y[ind]
        yerr_dn_all = abs(y[ind] - ydown[ind])
        yerr_up_all = abs(y[ind] - yup[ind])

        ind = np.where(binm == bin_mass)
        xin = x[ind]
        if(delta == False):
           yin = y[ind]
           yerr_dn = abs(yin - ydown[ind])
           yerr_up = abs(yin - yup[ind])
        else:
           nbins = len(binm[ind])
           yin = y[ind] - yin_all[0:nbins]
           yerr_dn = np.sqrt((y[ind] - ydown[ind])**2 + yerr_dn_all[0:nbins]**2)
           yerr_up = np.sqrt((y[ind] - yup[ind])**2 + yerr_up_all[0:nbins]**2)

        ax.errorbar(xin, yin, yerr=[yerr_dn, yerr_up], ls='None', mfc=color, fillstyle='full', ecolor = color, mec=color,marker='o',markersize=7, label="Brown+17" if lab else None)


    for i, sp in enumerate(subp):
        ax = fig.add_subplot(sp)
        plt.subplots_adjust(bottom=0.15, left=0.15)
        ytitle = "" 
        if(i == 0):
           ytitle = ytit
        prepare_ax(ax, 9, 11.5, -4, 1, xtit, ytitle)
       
        #Predicted SMHM
        ind = np.where(mh1_relation_satellites_halos[0,i,1,:] != 0)
        xplot = mh1_relation_satellites_halos[0,i,0,ind]
        yplot = mh1_relation_satellites_halos[0,i,1,ind]
        ax.plot(xplot[0],yplot[0],color='r', linestyle='solid', label=labels[i])
        ind = np.where(mh1_relation_satellites_halos[1,i,1,:] != 0)
        xplot = mh1_relation_satellites_halos[1,i,0,ind]
        yplot = mh1_relation_satellites_halos[1,i,1,ind]
        ax.plot(xplot[0],yplot[0],color='r', linestyle='dashed', label='with sat/cen confusion')

        if(i == 0):
           plot_obs_brown17(ax, bin_mass= i, mass=True, lab = True)
        else:
           plot_obs_brown17(ax, bin_mass= i, mass=True, lab = False)

        common.prepare_legend(ax, ['r','r','Salmon'])

    plt.tight_layout()
    common.savefig(output_dir, fig, "HI_stacking_satellites.pdf")


    fig = plt.figure(figsize=(16,4.5))
    xtit="$\\rm log_{10} (\\rm sSFR/yr^{-1})$"

    for i, sp in enumerate(subp):
        ax = fig.add_subplot(sp)
        plt.subplots_adjust(bottom=0.15, left=0.15)
 
        ytitle = ""        
        if(i == 0):
           ytitle = ytit
        prepare_ax(ax, -13, -8, -4, 1, xtit, ytitle)
       
        #Predicted SMHM
        ind = np.where(mh1_relation_ssfr_satellites_halos[0,i,1,:] != 0)
        xplot = mh1_relation_ssfr_satellites_halos[0,i,0,ind]
        yplot = mh1_relation_ssfr_satellites_halos[0,i,1,ind]
        ax.plot(xplot[0],yplot[0],color='r', linestyle='solid', label=labels[i])
        ind = np.where(mh1_relation_ssfr_satellites_halos[1,i,1,:] != 0)
        xplot = mh1_relation_ssfr_satellites_halos[1,i,0,ind]
        yplot = mh1_relation_ssfr_satellites_halos[1,i,1,ind]
        ax.plot(xplot[0],yplot[0],color='r', linestyle='dashed', label='with sat/cen confusion')

        if(i == 0):
           plot_obs_brown17(ax, bin_mass= i, mass=False, lab = True)
        else:
           plot_obs_brown17(ax, bin_mass= i, mass=False, lab = False)

        common.prepare_legend(ax, ['r','r','Salmon'])

    plt.tight_layout()
    common.savefig(output_dir, fig, "HI_stacking_satellites_ssfr.pdf")


    def plot_lagos18_HIstacking_relativediff(ax, mass_bin=1, mass=True, color='k'):
        if(mass == False):
            m, reldiff = common.load_observation(obs_dir, 'Models/SharkVariations/HIStackingSatellite_SSFR_RelativeDiff.dat', (0, 1))
            if(mass_bin == 1):
                i_i = 0
                i_f = 9 
            elif(mass_bin == 2):
                i_i = 9
                i_f = 19
            elif(mass_bin == 3):
                i_i = 19
                i_f = 29
            elif(mass_bin == 4):
                i_i=29
                i_f=39
        else:
            m, reldiff = common.load_observation(obs_dir, 'Models/SharkVariations/HIStackingSatellite_Mstar_RelativeDiff.dat', (0, 1))
            if(mass_bin == 1):
                i_i = 0
                i_f = 4 
            elif(mass_bin == 2):
                i_i = 4
                i_f = 9
            elif(mass_bin == 3):
                i_i = 9
                i_f = 15
            elif(mass_bin == 4):
                i_i=15
                i_f=21
        ax.plot(m[i_i:i_f], reldiff[i_i:i_f], linestyle='dashed', color=color)


    ###################################
    #   Plots global mass densities
    fig = plt.figure(figsize=(5,4.5))

    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm \\Delta\\,log_{10}\\langle M_{\\rm HI}/M_{\\star}\\rangle$"


    colors = ['k','blue', 'darkgreen', 'orange', 'red']

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)
    common.prepare_ax(ax, 9, 12, -1, 0.7, xtit, ytit, locators=(0.5, 0.5, 0.2, 0.2))

    for i, lab in enumerate(labels):
        if(i >= 1):
           plot_obs_brown17(ax, bin_mass= i, mass=True, lab = False, delta=True, color=colors[i])
           plot_lagos18_HIstacking_relativediff(ax, mass_bin=i, mass=True, color=colors[i])
           ind = np.where((mh1_relation_satellites_halos[0,i,1,:] != 0) & (np.isinf(mh1_relation_satellites_halos[0,i,1,:]) == False))
           xplot = mh1_relation_satellites_halos[0,i,0,ind] 
           yplot = mh1_relation_satellites_halos[0,i,1,ind] - mh1_relation_satellites_halos[0,0,1,ind]
           ax.plot(xplot[0],yplot[0],color=colors[i], linestyle='solid', label=labels[i])

           #ind = np.where((mh1_relation_satellites_halos[1,i,1,:] != 0) & (np.isinf(mh1_relation_satellites_halos[1,i,1,:]) == False))
           #xplot = mh1_relation_satellites_halos[1,i,0,ind]
           #yplot = mh1_relation_satellites_halos[1,i,1,ind] - mh1_relation_satellites_halos[1,0,1,ind]
           #ax.plot(xplot[0],yplot[0],color=colors[i], linestyle='dashed', label='with sat/cen confusion')
    common.prepare_legend(ax, ['blue', 'darkgreen', 'orange', 'red'], loc='best')

    plt.tight_layout()
    common.savefig(output_dir, fig, "HI_stacking_satellites_RelativeDiff.pdf")

    ###################################
    #   Plots global mass densities
    fig = plt.figure(figsize=(5,4.5))

    ytit="$\\rm \\Delta\\,log_{10}\\langle M_{\\rm HI}/M_{\\star}\\rangle$"
    xtit="$\\rm log_{10} (\\rm sSFR/yr^{-1})$"

    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.15, left=0.15)
    common.prepare_ax(ax, -12, -9, -1, 1, xtit, ytit, locators=(0.5, 0.5, 0.2, 0.2))

    for i, lab in enumerate(labels):
        if(i >= 1):
           plot_obs_brown17(ax, bin_mass= i, mass=False, lab = False, delta=True, color=colors[i])
           plot_lagos18_HIstacking_relativediff(ax, mass_bin=i, mass=False, color=colors[i])
           ind = np.where((mh1_relation_ssfr_satellites_halos[0,i,1,:] != 0) & (np.isinf(mh1_relation_ssfr_satellites_halos[0,i,1,:]) == False))
           xplot = mh1_relation_ssfr_satellites_halos[0,i,0,ind] 
           yplot = mh1_relation_ssfr_satellites_halos[0,i,1,ind] - mh1_relation_ssfr_satellites_halos[0,0,1,ind]
           ax.plot(xplot[0],yplot[0],color=colors[i], linestyle='solid')
           #ind = np.where((mh1_relation_satellites_halos[1,i,1,:] != 0) & (np.isinf(mh1_relation_satellites_halos[1,i,1,:]) == False))
           #xplot = mh1_relation_satellites_halos[1,i,0,ind]
           #yplot = mh1_relation_satellites_halos[1,i,1,ind] - mh1_relation_satellites_halos[1,0,1,ind]
           #ax.plot(xplot[0],yplot[0],color=colors[i], linestyle='dashed', label='with sat/cen confusion')
    #common.prepare_legend(ax, ['blue', 'darkgreen', 'orange', 'red'], loc='best')

    plt.tight_layout()
    common.savefig(output_dir, fig, "HI_stacking_satellites_ssfr_RelativeDiff.pdf")



def plot_molecular_gas_fraction(plt, output_dir, obs_dir, mgas_gals, mgas_relation, mh1_gals, mh1_relation, mh2_gals, mh2_relation, 
    mgas_relation_ltg, mh2_relation_ltg, mh1_relation_ltg, mgas_relation_etg, mh2_relation_etg, mh1_relation_etg, 
    mgas_ms_relation_ltg, mh2_ms_relation_ltg, mh1_ms_relation_ltg, mgas_ms_relation_etg, mh2_ms_relation_etg, 
    mh1_ms_relation_etg):

    xmin, xmax, ymin, ymax = 9, 12, -3, 1
    fig = plt.figure(figsize=(11,11))

    # First subplot
    ax = fig.add_subplot(321)
    plt.subplots_adjust(left=0.15)

    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm HI+H_2}/M_{\\star})$"
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    #Predicted relation for all galaxies
    ind = np.where((mgas_gals[0,:] > 0) & (mgas_gals[1,:] != 0) )
    xdata = mgas_gals[0,ind]
    ydata = mgas_gals[1,ind]
    us.density_contour(ax, xdata[0], ydata[0], 30, 30) #, **contour_kwargs)

    def plot_mrelation(mrelation, color, label=None, linestyle=None):
        ind = np.where(mrelation[0,:] != 0)
        xplot = xmf[ind]
        yplot = mrelation[0,ind]
        linestyle = linestyle or ''
        ax.plot(xplot,yplot[0], color=color, label=label, linestyle=linestyle)

    def plot_mrelation_fill(mrelation, color, colorfill, label=None, linestyle=None):
        ind = np.where(mrelation[0,:] != 0)
        xplot = xmf[ind]
        yplot = mrelation[0,ind]
        errdn = mrelation[1,ind]
        errup = mrelation[2,ind]

        ax.plot(xplot,yplot[0], color=color, label=label, linestyle=linestyle)
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor=colorfill, alpha=0.2,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor=colorfill, alpha=0.2,interpolate=True)

    plot_mrelation(mgas_relation, 'k', linestyle='solid', label="Shark all galaxies")

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    add_observations_to_plot(obs_dir, 'NeutralGasRatio_NonDetEQZero.dat', ax, '^', "xCOLDGAS+xGASS(0)", color='grey')
    add_observations_to_plot(obs_dir, 'NeutralGasRatio_NonDetEQUpperLimits.dat', ax, 'v', "xCOLDGAS+xGASS(UL)")

    common.prepare_legend(ax, ['k','k','k'])

    # Second subplot
    ax = fig.add_subplot(322)
    xmin, xmax, ymin, ymax = 9, 12, -4.5, 1

    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    plot_mrelation_fill(mgas_relation_ltg, 'b', 'b',label="Shark LTGs $(\\rm B/T)_{\\rm bar}$",linestyle='solid')
    plot_mrelation_fill(mgas_relation_etg, 'r', 'r',label="Shark ETGs $(\\rm B/T)_{\\rm bar}$",linestyle='solid')
    plot_mrelation(mgas_ms_relation_ltg, 'b',label="Shark LTGs $(\\rm B/T)_{\star}$",linestyle='dotted')
    plot_mrelation(mgas_ms_relation_etg, 'r',label="Shark ETGs $(\\rm B/T)_{\star}$",linestyle='dotted')

    # Legend
    common.prepare_legend(ax, ['b','r','b','r','k'],loc=3)

    # Third subplot
    ax = fig.add_subplot(323)
    plt.subplots_adjust(left=0.15)

    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm HI}/M_{\\star})$"
    xmin, xmax, ymin, ymax = 9, 12, -3, 1
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    #Predicted relation
    ind = np.where((mh1_gals[0,:] > 0) & (mh1_gals[1,:] != 0) )
    xdata = mh1_gals[0,ind]
    ydata = mh1_gals[1,ind]
    us.density_contour(ax, xdata[0], ydata[0], 30, 30) #, **contour_kwargs)
    plot_mrelation(mh1_relation, 'k', linestyle='solid')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    add_observations_to_plot(obs_dir, 'HIGasRatio_NonDetEQZero.dat', ax, '^', "xGASS(0)", color='grey')
    add_observations_to_plot(obs_dir, 'HIGasRatio_NonDetEQUpperLimits.dat', ax, 'v', "xGASS(UL)")

    x, y, yerr_down, yerr_up = common.load_observation(obs_dir, 'Gas/Parkash18.dat', (0, 1, 2, 3))
    ax.errorbar(x,y-x,yerr=[(y-x) - (yerr_down-x),(yerr_up-x) - (y-x)], ls='None', mfc='r', fillstyle='full', ecolor = 'r', mec='r',marker='s',markersize=7, label="Parkash+18")

    m, mrat, merr = common.load_observation(obs_dir, 'Gas/RHI-Mstars_Brown15.dat', [0,1,2])
    errdn = np.log10(mrat) - np.log10(mrat - merr) 
    errup = np.log10(mrat + merr) - np.log10(mrat) 
    ax.errorbar(m,np.log10(mrat),yerr=[errdn,errup], ls='None', mfc='Salmon', fillstyle='full', ecolor = 'Salmon', mec='Salmon',marker='o',markersize=7, label="Brown+15")

    # Legend
    common.prepare_legend(ax, ['k','k','r','Salmon'], loc=1)

    # Fourth subplot
    ax = fig.add_subplot(324)
    xmin, xmax, ymin, ymax = 9, 12, -4.5, 1
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    plot_mrelation_fill(mh1_relation_ltg, 'b', 'b',linestyle='solid')
    plot_mrelation_fill(mh1_relation_etg, 'r', 'r',linestyle='solid')
    plot_mrelation(mh1_ms_relation_ltg,  'b',linestyle='dotted')
    plot_mrelation(mh1_ms_relation_etg,  'r',linestyle='dotted')

    add_observations_to_plot(obs_dir, 'RHI-Mstars_Callette18-LTGs.dat', ax, 's', "Calette+18 LTGs", color='grey', err_absolute=True)
    add_observations_to_plot(obs_dir, 'RHI-Mstars_Callette18-ETGs.dat', ax, 'o', "Calette+18 ETGs", color='grey', err_absolute=True)

    # Legend
    common.prepare_legend(ax, ['grey','grey','grey'],loc=1)

    # Fifth subplot
    ax = fig.add_subplot(325)
    plt.subplots_adjust(left=0.15)

    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm H_2}/M_{\\star})$"
    xmin, xmax, ymin, ymax = 9, 12, -3, 1
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    #Predicted relation
    ind = np.where((mh2_gals[0,:] > 0) & (mh2_gals[1,:] != 0) )
    xdata = mh2_gals[0,ind]
    ydata = mh2_gals[1,ind]
    us.density_contour(ax, xdata[0], ydata[0], 30, 30) #, **contour_kwargs)
    plot_mrelation(mh2_relation, 'k', linestyle='solid')

    #Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
    add_observations_to_plot(obs_dir, 'MolecularGasRatio_NonDetEQZero.dat', ax, '^', "xCOLDGASS(0)", color='grey')
    add_observations_to_plot(obs_dir, 'MolecularGasRatio_NonDetEQUpperLimits.dat', ax, 'v', "xCOLDGASS(UL)")

    common.prepare_legend(ax, ['k','k','k'], loc = 1)

    # Fourth subplot
    ax = fig.add_subplot(326)
    xmin, xmax, ymin, ymax = 9, 12, -4.5, 1
    prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit)

    plot_mrelation_fill(mh2_relation_ltg, 'b', 'b',linestyle='solid')
    plot_mrelation_fill(mh2_relation_etg, 'r', 'r',linestyle='solid')
    plot_mrelation(mh2_ms_relation_ltg, 'b',linestyle='dotted')
    plot_mrelation(mh2_ms_relation_etg, 'r',linestyle='dotted')

    add_observations_to_plot(obs_dir, 'RH2-Mstars_Callette18-LTGs.dat', ax, 's', "Calette+18 LTGs",color='grey', err_absolute=True)
    add_observations_to_plot(obs_dir, 'RH2-Mstars_Callette18-ETGs.dat', ax, 'o', "Calette+18 ETGs",color='grey', err_absolute=True)

    # Legend
    common.prepare_legend(ax, ['grey','grey','grey'],loc=1)

    common.savefig(output_dir, fig, "molecular_gas_fraction.pdf")

def plot_h1h2_gas_fraction(plt, output_dir, mhr_relation, mhr_relation_cen, mhr_relation_sat):

    fig = plt.figure(figsize=(5,4.5))

    ax = fig.add_subplot(111)
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10}(M_{\\rm H_2}/M_{\\rm HI})$"
    prepare_ax(ax, 8, 12, -3, 1.0, xtit, ytit)

    # Predicted SMHM
    ind = np.where(mhr_relation[0,:] != 0)
    xplot = xmf[ind]
    yplot = mhr_relation[0,ind]
    errdn = mhr_relation[1,ind]
    errup = mhr_relation[2,ind]

    ax.errorbar(xplot,yplot[0],color='k', label="all galaxies")
    ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'k', mec='k',marker='+',markersize=2)


    ind = np.where(mhr_relation_cen[0,:] != 0)
    xplot = xmf[ind]
    yplot = mhr_relation_cen[0,ind]
    ax.errorbar(xplot,yplot[0],color='b',linestyle='dotted', label="centrals")

    ind = np.where(mhr_relation_sat[0,:] != 0)
    xplot = xmf[ind]
    yplot = mhr_relation_sat[0,ind]
    ax.errorbar(xplot,yplot[0],color='r',linestyle='dashed', label="satelites")

    common.prepare_legend(ax, ['k','b','r','grey','grey'])
    common.savefig(output_dir, fig, "HIH2_gas_fraction.pdf")


def main(model_dir, output_dir, redshift_table, subvols, obs_dir):


    zlist = [0, 0.381963715160695]

    plt = common.load_matplotlib()
    fields = {'galaxies': ('type', 'mstars_disk', 'mstars_bulge',
                           'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge', 'mvir_hosthalo', 
                           'sfr_disk', 'sfr_burst','id_halo_tree')}

    for index, snapshot in enumerate(redshift_table[zlist]):

         hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)

         (mgas_relation, mgas_relation_cen, mgas_relation_sat,
         mh2_gals, mh1_gals, mgas_gals,
         mh2_relation, mh1_relation, mhr_relation, mhr_relation_cen, mhr_relation_sat,
         mgas_relation_ltg, mh2_relation_ltg, mh1_relation_ltg,
         mgas_relation_etg, mh2_relation_etg, mh1_relation_etg,
         mgas_ms_relation_ltg, mh2_ms_relation_ltg, mh1_ms_relation_ltg,
         mgas_ms_relation_etg, mh2_ms_relation_etg, mh1_ms_relation_etg, 
         mh1_relation_satellites_halos, mh1_relation_ssfr_satellites_halos,
         mh1_relation_satellites_halos_stack, mh1_relation_ssfr_satellites_halos_stack) = prepare_data(index, zlist[index], hdf5_data)

         if(index == 0):
            plot_cold_gas_fraction(plt, output_dir, obs_dir, mgas_relation, mgas_relation_cen, mgas_relation_sat)

            plot_HI_stacking(plt, output_dir, obs_dir, mh1_relation_satellites_halos_stack, mh1_relation_ssfr_satellites_halos_stack)
            
            plot_molecular_gas_fraction(plt, output_dir, obs_dir, mgas_gals, mgas_relation, mh1_gals, mh1_relation, mh2_gals, mh2_relation, mgas_relation_ltg, 
                mh2_relation_ltg, mh1_relation_ltg, mgas_relation_etg, mh2_relation_etg, mh1_relation_etg, mgas_ms_relation_ltg, mh2_ms_relation_ltg, 
                mh1_ms_relation_ltg, mgas_ms_relation_etg, mh2_ms_relation_etg, mh1_ms_relation_etg)
            
            plot_h1h2_gas_fraction(plt, output_dir, mhr_relation, mhr_relation_cen, mhr_relation_sat)

if __name__ == '__main__':
    main(*common.parse_args())
