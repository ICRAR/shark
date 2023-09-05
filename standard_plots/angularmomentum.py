
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
"""Angular momentum plots"""

import functools

import numpy as np

import common
import utilities_statistics as us

# Initialize arguments
zlist = [0, 0] #, 1, 2)

##################################
#Constants
RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun
XH = 1.33

#stellar mass bins
mlow = 6.35
mupp = 12.5
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0
#stellar mass bins for observations
mlowo = 6.45
muppo = 12.5
dmo = 0.22
mbinso = np.arange(mlowo,muppo,dmo)
xmfo = mbinso + dmo/2.0


#DM bins
mlowh = 10.0
mupph = 15.0
mbinsh = np.arange(mlowh,mupph,dm)
xmfh   = mbinsh + dm/2.0

dmobs = 0.4
mbins_obs = np.arange(mlow,mupp,dmobs)
xmf_obs = mbins_obs + dmobs/2.0

llow = -3.0
lupp = 6.0
dl   = 0.15
lbins = np.arange(llow,lupp,dl)
xlf   = lbins + dl/2.0


def prepare_data(hdf5_data, index, sam_stars_disk, sam_gas_disk_atom, sam_gas_disk_atom2, sam_gas_disk_mol, sam_gas_disk_atom_ms, sam_halo, sam_ratio_halo_disk, sam_ratio_halo_gal, 
                 sam_ratio_halo_disk_gas, disk_size_sat, disk_size_cen, bulge_size, sam_vs_sam_halo_disk, sam_vs_sam_halo_gal,
                 sam_vs_sam_halo_disk_gas, sam_bar, sam_stars, sam_stars2, vmax_halo_gal, sam_barv2, disk_size, rcomb): #, phot_data, phot_data_nod, nbands):


    bin_it = functools.partial(us.wmedians, xbins=xmf, low_numbers=False)
    bin_it_largem = functools.partial(us.wmedians, xbins=xmf_obs, low_numbers=False)
    bin_it_halo = functools.partial(us.wmedians, xbins=xmfh, low_numbers=True)
    bin_it_j    = functools.partial(us.wmedians, xbins=xlf, low_numbers=True)

    (h0, _, mdisk, mbulge, mburst_mergers, mburst_diskins, mstars_bulge_mergers_assembly, mstars_bulge_diskins_assembly, 
     mBH, rdisk, rbulge, rg_disk, rg_bulge, typeg, specific_angular_momentum_disk_star, specific_angular_momentum_bulge_star, 
     specific_angular_momentum_disk_gas, specific_angular_momentum_bulge_gas, specific_angular_momentum_disk_gas_atom, 
     specific_angular_momentum_disk_gas_mol, lambda_sub, mvir_s, mvir, matom_disk, mmol_disk, mgas_disk,
     matom_bulge, mmol_bulge, mgas_bulge, sfr_disk, sfr_bulge, vmax, lx, ly, lz) = hdf5_data

    ind = np.where(mmol_disk > 0) 
    print("sAM mol when mgas_mol>0", specific_angular_momentum_disk_gas_mol[ind])
    print("number of galaxies with mmol_disk > 0", len(mmol_disk[ind]))
    ind = np.where((mmol_disk > 0) & (specific_angular_momentum_disk_gas_mol == 0) & (typeg ==0))
    print("number of galaxies with mmol_disk > 0 but jmol_disk == 0", len(mmol_disk[ind]), np.median(mmol_disk[ind]), np.median(mdisk[ind]+mbulge[ind]))
    ind = np.where((mmol_disk > 0) & (specific_angular_momentum_disk_gas_mol > 0) & (typeg ==0)) 
    print("number of galaxies with mmol_disk > 0 and jmol_disk > 0", len(mmol_disk[ind]), np.median(mmol_disk[ind]), np.median(mdisk[ind]+mbulge[ind]))


    sam_subhalo_fromL = np.sqrt(lx**2.0 + ly**2.0 + lz**2.0) * 1e3 / mvir_s #in km/s * kpc
    lambda_fromL = np.sqrt(lx**2.0 + ly**2.0 + lz**2.0) / mvir_s * 1.5234153 / (4.3e-9 * mvir_s)**0.666 * 67.77**0.33;

    #ind = np.where((mvir_s > 0) & (typeg==0))
    #print(max(lambda_fromL[ind]), np.median(lambda_fromL[ind]))

    sam_subhalo_rel = bin_it_halo(x=np.log10(mvir_s[ind]/h0),
                             y=np.log10(sam_subhalo_fromL[ind]/h0)) #specific_angular_momentum_disk_star[ind]) - np.log10(float(h0)))


    ind = np.where(mdisk + mbulge > 0)
    #SEDs_dust = np.zeros(shape = (len(mdisk[ind]), 5, nbands))
    #SEDs_nodust = np.zeros(shape = (len(mdisk[ind]), 5, nbands))
    mstartot = mdisk[ind] + mbulge[ind]
    rgal_star = (rdisk * mdisk + rbulge * mbulge ) / ( mdisk + mbulge)
    rgal_gas  = (rg_disk * mgas_disk + rg_bulge * mgas_bulge ) / ( mgas_disk + mgas_bulge)
   
    ind = np.where(mdisk+mbulge > 0)
    rcomb[index,:] = bin_it_largem(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                            y=np.log10((mdisk[ind]*rdisk[ind]  + mbulge[ind]*rbulge[ind])*MpcToKpc / (mdisk[ind]+mbulge[ind]))- np.log10(float(h0)))


    mbulge_mergers = mburst_mergers + mstars_bulge_mergers_assembly
    zero_bulge = np.where(rbulge <= 0)
    if(len(rbulge) == len(rbulge[zero_bulge])):
            #case where there is zero bulge build up.
            rbulge[zero_bulge] = 1e-10
            specific_angular_momentum_bulge_star[zero_bulge] = 1.0
            mbulge[zero_bulge] = 10.0


    H           = h0*100.0 * pow(0.3121*pow(1+zlist[index],3.0) + 0.6879, 0.5)
    sam_subhalo = 1.41421356237 * G**0.666 * lambda_sub * mvir_s**0.666 / (H)**0.333
    sam_hhalo   = 1.41421356237 * G**0.666 * lambda_sub * mvir**0.666 / (H)**0.333

    #calculate effective specific angular momentum by mass-weighting the contributions from the disk and bulge
    jstars      = (specific_angular_momentum_disk_star * mdisk + specific_angular_momentum_bulge_star * mbulge) / (mdisk+mbulge)
    jbar        = (specific_angular_momentum_disk_star * mdisk + specific_angular_momentum_disk_gas * mgas_disk ) / (mdisk + mbulge + mgas_disk + mgas_bulge)  
                   #specific_angular_momentum_bulge_star * mbulge + specific_angular_momentum_bulge_gas * mgas_bulge) / (mdisk + mbulge + mgas_disk + mgas_bulge)
    mbar        = (mdisk + mbulge + mgas_disk + mgas_bulge)
    jbarv2      = (specific_angular_momentum_disk_star * mdisk + specific_angular_momentum_disk_gas * mgas_disk + specific_angular_momentum_bulge_star * mbulge + specific_angular_momentum_bulge_gas * mgas_bulge) / (mdisk + mbulge + mgas_disk + mgas_bulge)
    #(specific_angular_momentum_disk_gas_atom * matom_disk + specific_angular_momentum_disk_star * mdisk + specific_angular_momentum_bulge_star * mbulge + matom_bulge * specific_angular_momentum_bulge_gas) / (mdisk + mbulge + matom_disk + matom_bulge)
    mbarv2      = (mdisk + mbulge + mgas_disk + mgas_bulge)
    jatom       = (specific_angular_momentum_disk_gas_atom * matom_disk + matom_bulge * specific_angular_momentum_bulge_gas) / (matom_disk + matom_bulge)
    jmol        = (specific_angular_momentum_disk_gas_mol * mmol_disk + mmol_bulge * specific_angular_momentum_bulge_gas) / (mmol_disk + mmol_bulge)

    jgas        = (specific_angular_momentum_disk_gas * mgas_disk + specific_angular_momentum_bulge_gas * mgas_bulge) / (mgas_bulge + mgas_disk)
    mgas        = (mgas_bulge + mgas_disk)
    vdisk  = specific_angular_momentum_disk_star / rdisk / 2.0 #in km/s
    vbulge = specific_angular_momentum_bulge_star / rbulge / 2.0 #in km/s

    specific_angular_momentum_disk = (specific_angular_momentum_disk_star * mdisk + specific_angular_momentum_disk_gas * mgas_disk) / (mdisk + mgas_disk)

    vr_halo = G**0.66 * mvir**0.66 / (H)**0.33
    lh = lambda_sub
    lj = np.zeros(shape = (4, len(mdisk))) 
    lm = np.zeros(shape = (4, len(mdisk))) 
    
    lj[0,:] = specific_angular_momentum_disk  / 1.41421356237 / vr_halo
    lj[1,:] = specific_angular_momentum_disk_star  / 1.41421356237 / vr_halo
    lj[2,:] = specific_angular_momentum_disk_gas   / 1.41421356237 / vr_halo
    lj[3,:] = jstars / 1.41421356237 / vr_halo

    lm[0,:] = specific_angular_momentum_disk  / 1.41421356237 / G**0.66 * (H)**0.33 / (mdisk + mgas_disk)**0.66
    lm[1,:] = specific_angular_momentum_disk_star  / 1.41421356237 / G**0.66 * (H)**0.33 / (mdisk)**0.66
    lm[2,:] = specific_angular_momentum_disk_gas   / 1.41421356237 / G**0.66 * (H)**0.33 / (mgas_disk)**0.66
    lm[3,:] = jstars  / 1.41421356237 / G**0.66 * (h0*100.0)**0.33 / (mdisk + mbulge)**0.66

    bt = np.zeros(shape = (len(mdisk)))
    ms = np.zeros(shape = (len(mdisk)))
    ssfr = np.zeros(shape = (len(mdisk)))
    ind = np.where(mdisk+mbulge > 0) 
    bt[ind] = mbulge[ind] / (mdisk[ind] + mbulge[ind])
    ms[ind] = np.log10(mdisk[ind] + mbulge[ind])
    ind = np.where((mdisk+mbulge > 0) & (sfr_disk + sfr_bulge > 0))
    ssfr[ind] = np.log10((sfr_disk[ind] + sfr_bulge[ind])) - ms[ind] #in Gyr

    thresh = [0, 0.25, 0.5, 0.75]
    thresh2 = [0, 0.25, 0.5, 0.75]

    mass_cut = 8.5
    for c in range(0,len(thresh)):
        ind = np.where((jstars > 0) & (ms > mass_cut) & (mdisk+mbulge > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c]))
        sam_stars[index,:,:,c]                = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                                  y=np.log10(jstars[ind]) - np.log10(float(h0))) #specific_angular_momentum_disk_star[ind]) - np.log10(float(h0)))
        ind = np.where((jstars > 0) & (ms > mass_cut)  & (mdisk+mbulge > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh2[c]))
        sam_stars2[index,:,:,c]               = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                                  y=np.log10(specific_angular_momentum_disk_star[ind]) - np.log10(float(h0)))
        ind = np.where((jatom > 0) & (ms > mass_cut)  & (mdisk+mbulge > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh2[c]))
        sam_gas_disk_atom_ms[index,:,:,c]     = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                                  y=np.log10(jatom[ind]) - np.log10(float(h0)))
        ind = np.where((jmol > 0) & (ms > mass_cut)  & (mdisk+mbulge > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh2[c]))
        sam_gas_disk_mol[index,:,:,c]         = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                                  y=np.log10(jmol[ind]) - np.log10(float(h0)))

        ind = np.where((specific_angular_momentum_disk_star > 0) & (mdisk+mbulge > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c]))
        sam_stars_disk[index,:,:,c]           = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
                                                  y=np.log10(specific_angular_momentum_disk_star[ind]) - np.log10(float(h0))) #specific_angular_momentum_disk_star[ind]) - np.log10(float(h0)))

        ind = np.where((jbar > 0) & (ms > mass_cut) & (mbar > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c]))
        sam_bar[index,:,:,c]                  = bin_it(x=np.log10(mbar[ind]) - np.log10(float(h0)),
                                                  y=np.log10(jbar[ind]) - np.log10(float(h0))) #specific_angular_momentum_disk_star[ind]) - np.log10(float(h0)))

        ind = np.where((jbarv2 > 0) & (ms > mass_cut) & (mbarv2 > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c])) 
        sam_barv2[index,:,:,c]                = bin_it(x=np.log10(mbarv2[ind]) - np.log10(float(h0)),
                                                  y=np.log10(jbarv2[ind]) - np.log10(float(h0))) #specific_angular_momentum_disk_star[ind]) - np.log10(float(h0)))

        ind = np.where((jstars > 0) & (ms > mass_cut)  & (mdisk+mbulge > 0) & (typeg >= 0) & (mdisk/(mdisk+mbulge) >= thresh[c]))
        sam_ratio_halo_gal[index,:,:,c]       = bin_it_halo(x=np.log10(mvir_s[ind]) - np.log10(float(h0)),
                                                  y=np.log10(jstars[ind]/sam_subhalo[ind]))
        sam_vs_sam_halo_gal[index,:,:,c]      = bin_it_j(x=np.log10(sam_subhalo[ind]),
                                                  y=np.log10(jstars[ind])  - np.log10(float(h0)))
        ind = np.where(typeg == 0)
        vmax_halo_gal[index,:,:,c]            = bin_it_halo(x=np.log10(mvir_s[ind]) - np.log10(float(h0)),
                                                  y=np.log10(vmax[ind]))

        ind = np.where((specific_angular_momentum_disk_star > 0) & (ms > mass_cut)  & (mdisk+mbulge > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c]))
        sam_ratio_halo_disk[index,:,:,c]      = bin_it_halo(x=np.log10(mvir_s[ind]) - np.log10(float(h0)),
                                                  y=np.log10(specific_angular_momentum_disk_star[ind]/sam_subhalo[ind]))
        sam_vs_sam_halo_disk[index,:,:,c]     = bin_it_j(x=np.log10(sam_subhalo[ind]),
                                                  y=np.log10(specific_angular_momentum_disk_star[ind])  - np.log10(float(h0)))
    
        ind = np.where((specific_angular_momentum_disk_gas > 0) & (ms > mass_cut)  & (mdisk+mbulge > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c]))
        sam_ratio_halo_disk_gas[index,:,:,c]  = bin_it_halo(x=np.log10(mvir_s[ind]) - np.log10(float(h0)),
                                                  y=np.log10(specific_angular_momentum_disk_gas[ind]/sam_subhalo[ind]))
        sam_vs_sam_halo_disk_gas[index,:,:,c] = bin_it_j(x=np.log10(sam_subhalo[ind]),
                                                  y=np.log10(specific_angular_momentum_disk_gas[ind])  - np.log10(float(h0)))
     
        ind = np.where((jatom > 0) & (ms > mass_cut) & (matom_disk+matom_bulge > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c]) & (ms < 12))
        sam_gas_disk_atom[index,:,:,c]        = bin_it(x=np.log10(matom_disk[ind] + matom_bulge[ind]) - np.log10(float(h0)),
                                                  y=np.log10(jatom[ind]) - np.log10(float(h0)))
        ind = np.where((jgas > 0) & (ms > mass_cut) & (mgas > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c]) & (ms < 12))
        sam_gas_disk_atom2[index,:,:,c]        = bin_it(x=np.log10(mgas[ind]) - np.log10(float(h0)),
                                                  y=np.log10(jgas[ind]) - np.log10(float(h0)))
      
        #ind = np.where((jmol > 0) & (ms > mass_cut)  & (mdisk+mbulge > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh2[c]))
        #sam_gas_disk_mol[index,:,:,c]         = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)),
        #                                          y=np.log10(specific_angular_momentum_disk_gas_mol[ind]) - np.log10(float(h0)))
      
        ind = np.where((sam_subhalo > 0) & (mdisk+mbulge > 0) & (ms > mass_cut)  & (typeg == 0) & (mdisk/(mdisk+mbulge) >= thresh[c]))
        sam_halo[index,:,:,c]                 = bin_it(x=np.log10(mdisk[ind]+mbulge[ind]) - np.log10(float(h0)), 
            		                          y=np.log10(sam_subhalo[ind]))
    

    ind = np.where((mdisk > 0) & (typeg == 0) & (mdisk/(mdisk+mbulge) > 0.5))
    disk_size_cen[index,:]  = bin_it(x=np.log10(mdisk[ind]) - np.log10(float(h0)),
                                    y=np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)))

    ind = np.where((mdisk > 0) & (typeg > 0) & (mdisk/(mdisk+mbulge) > 0.5))
    disk_size_sat[index,:] = bin_it(x=np.log10(mdisk[ind]) - np.log10(float(h0)),
                                    y=np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)))

    ind = np.where((mbulge > 0) & (mbulge/(mbulge+mdisk) > 0.5) & (rbulge > 1e-6))
    bulge_size[index,:] = bin_it_largem(x=np.log10(mbulge[ind]) - np.log10(float(h0)),
                                 y=np.log10(rbulge[ind]*MpcToKpc) - np.log10(float(h0)))

    ind = np.where((mbulge > 0) & (mdisk/(mbulge+mdisk) <= 0.5) & (rbulge > 1e-6))
    disk_size[index,:] = bin_it_largem(x=np.log10(mdisk[ind]) - np.log10(float(h0)),
                                 y=np.log10(rdisk[ind]*MpcToKpc) - np.log10(float(h0)))

    return (lh, lj, lm, bt, ms, ssfr, mass_cut)


def plot_sizes(plt, outdir, obsdir, disk_size_cen, disk_size_sat, bulge_size, vmax_halo_gal, disk_size, rcomb):

    #print ('sizes disks')
    #for i,j,p in zip(disk_size[0,0,:],disk_size[0,1,:],disk_size[0,2,:]):
    #    print( i,j,p)
    #print('sizes bulges')
    #for i,j,p in zip(bulge_size[0,0,:],bulge_size[0,1,:],bulge_size[0,2,:]):
    #    print(i,j,p)

    def plot_profuse_sizes(ax, component=0, color='DarkTurquoise', include_ambig = True, r_bulge_assumed = 1):
        bin_it = functools.partial(us.wmedians, xbins=xmfo, low_numbers=False)
        sizes_table = np.loadtxt(obsdir + '/SizeMass/' + 'ProFuse_cat_BestModel_v1.0.csv', delimiter=',')
        r=0
        mass=0
        if((component ==0) & (include_ambig == True)):
           #select disks of interest:
           ind = np.where((sizes_table[:,3] == 1) | (sizes_table[:,2] == 1) | (sizes_table[:,2] == 2) | (sizes_table[:,3] ==3)) #Disk or D+B or P+D or Ambg
           mass = sizes_table[ind,9]
           r = sizes_table[ind,22]
           lab = '(disks + ambig)'
           fill = 'full'
        elif((component ==0) & (include_ambig == False)):
           #select disks of interest:
           ind = np.where((sizes_table[:,3] == 1) | (sizes_table[:,2] == 1) | (sizes_table[:,2] == 2)) #Disk or D+B or P+D
           mass = sizes_table[ind,9]
           r = sizes_table[ind,22]
           lab = '(disks)'
           fill = 'none'
        elif((component ==1) & (include_ambig == True)):
           #select bulges of interest:
           ind = np.where((sizes_table[:,3] == 2) | (sizes_table[:,2] == 1) | (sizes_table[:,3] ==3) | (sizes_table[:,2] == 3)) #Sph or D+B or Ambg or P+D
           mass = sizes_table[ind,32]
           r = sizes_table[ind,45]
           ind = np.where(r == -99)
           r[ind] = r_bulge_assumed
           lab = '(b+Sph+ambig)'
           fill = 'full'
        elif((component ==1) & (include_ambig == False)):
           ind = np.where((sizes_table[:,3] == 2) | (sizes_table[:,2] == 1) | (sizes_table[:,2] == 3)) #Sph or D+B or P+D
           mass = sizes_table[ind,32]
           r = sizes_table[ind,45]
           ind = np.where(r == -99)
           r[ind] = r_bulge_assumed
           lab = '(b+Sph)'
           fill = 'none'
        elif(component == 4):
           #select bulges of interest:
           ind = np.where((sizes_table[:,3] == 2) | (sizes_table[:,3] ==3)) #Sph
           mass = sizes_table[ind,32]
           r = sizes_table[ind,45]
           ind = np.where(r == -99)
           r[ind] = r_bulge_assumed
           lab = '(Sph+ambig)'
           fill = 'none'
        elif(component == 5):
           #select bulges of interest:
           ind = np.where((sizes_table[:,2] == 1)) #Sph
           mass = sizes_table[ind,32]
           r = sizes_table[ind,45]
           ind = np.where(r == -99)
           r[ind] = r_bulge_assumed
           lab = '(b)'
           fill = 'none'

        elif(component == 2):
           mass = np.zeros(shape = len(sizes_table[:,3]))
           r = np.zeros(shape = len(sizes_table[:,3]))
           ind = np.where(sizes_table[:,3] == 2) #pure spheroids
           mass[ind] = sizes_table[ind,32]
           r[ind] = sizes_table[ind,45]
           ind = np.where(sizes_table[:,3] == 1) #pure disks
           mass[ind] = sizes_table[ind,9]
           r[ind] = sizes_table[ind,22]
           ind = np.where(sizes_table[:,3] == 3) #ambigious
           mass[ind] = sizes_table[ind,9]
           r[ind] = sizes_table[ind,22]
           ind = np.where(sizes_table[:,2] == 1) #B+D
           mass[ind] = np.log10(10**sizes_table[ind,32] + 10**sizes_table[ind,9])
           r[ind] = (sizes_table[ind,22]* 10**sizes_table[ind,9] + sizes_table[ind,45] * 10**sizes_table[ind,32]) / (10**sizes_table[ind,32] + 10**sizes_table[ind,9])
           ind = np.where(sizes_table[:,2] == 3) #P+D (here assume bulge = r_bulge_assumed [kpc])
           mass[ind] = np.log10(10**sizes_table[ind,32] + 10**sizes_table[ind,9])
           r[ind] = (sizes_table[ind,22]* 10**sizes_table[ind,9] + r_bulge_assumed * 10**sizes_table[ind,32]) / (10**sizes_table[ind,32] + 10**sizes_table[ind,9])
           lab = '(mass-weigh)'
           fill = 'full'

        ind = np.where(r > 0)
        mzr = bin_it(x=mass[ind], y=np.log10(r[ind]))
        ind = np.where(mzr[0,:] != 0)
        if(len(xmf[ind]) > 0):
            xplot = xmfo[ind]
            yplot = mzr[0,ind]
            errdn = mzr[1,ind]
            errup = mzr[2,ind]
            ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], color=color,linestyle='None', marker='s', fillstyle=fill, label='Bellstedt+23 ' + lab)
 

    rb, r16, r84 = common.load_observation(obsdir, 'Models/SharkVariations/SizeDisksAndBulges_OtherModels.dat', [0,1,2])
    
    fig = plt.figure(figsize=(5.5,12))
    xtit = "$\\rm log_{10} (\\rm M_{\\star,disk}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm r_{\\star,disk}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.3, 2
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # LTG ##################################
    ax = fig.add_subplot(311)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(10.6,0.1,'galaxy disks',fontsize=12)

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(disk_size[0,0,:] != 0)
    xplot = xmf_obs[ind]
    yplot = disk_size[0,0,ind]
    errdn = disk_size[0,1,ind]
    errup = disk_size[0,2,ind]
    ax.plot(xplot,yplot[0],color='b',linestyle='solid',label="Shark v2.0")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.3,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.3,interpolate=True)

    rdisk_am   = rb[30:59]
    rdisk_am16 = r16[30:59]
    rdisk_am84 = r84[30:59]
    ind = np.where(rdisk_am != 0)
    xplot = xmf[ind]
    yplot = rdisk_am[ind]
    errdn = rdisk_am16[ind]
    errup = rdisk_am84[ind]
    ax.plot(xplot,yplot,color='DodgerBlue',linestyle='dashed', label="Shark v1.1 (L18)")
    ax.fill_between(xplot,yplot,yplot-errdn, facecolor='DeepSkyBlue', linestyle='solid', alpha=0.4,interpolate=True)
    ax.fill_between(xplot,yplot,yplot+errup, facecolor='DeepSkyBlue', linestyle='solid', alpha=0.4,interpolate=True)

    #Lange et al. (2016)
    m,r = common.load_observation(obsdir, 'SizesAndAM/rdisk_L16.dat', [0,1])
    m = np.log10(m)
    r = np.log10(r)
    ax.plot(m[0:39], r[0:39], linestyle='dotted',color='k')
    ax.plot(m[39:78], r[39:78], linestyle='dotted',color='k')
    ax.plot(m[78:len(r)], r[78:len(r)], linestyle='dotted',color='k',label='Lange+16 (disks)')
    plot_profuse_sizes(ax, component=0, color='CadetBlue', include_ambig = True)
    plot_profuse_sizes(ax, component=0, color='darkBlue', include_ambig = False)

    common.prepare_legend(ax, ['b','DeepSkyBlue','black', 'CadetBlue','darkblue'], loc=2)

    # ETGs ##################################
    xtit = "$\\rm log_{10} (\\rm M_{\\star,bulge}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm r_{\\star,bulge}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.5, 2
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(312)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(10.6,-0.2,'galaxy bulges',fontsize=12)

    #Predicted size-mass for bulges in bulge-dominated systems
    ind = np.where(bulge_size[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf_obs[ind]
        yplot = bulge_size[0,0,ind]
        errdn = bulge_size[0,1,ind]
        errup = bulge_size[0,2,ind]
        ax.plot(xplot,yplot[0],color='r',linestyle='solid', label='Shark v2.0')
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.3,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.3,interpolate=True)

    rdisk_am   = rb[60:99]
    rdisk_am16 = r16[60:99]
    rdisk_am84 = r84[60:99]

    ind = np.where(rdisk_am != 0)
    xplot = xmf[ind]
    yplot = rdisk_am[ind]
    errdn = rdisk_am16[ind]
    errup = rdisk_am84[ind]
    ax.plot(xplot,yplot,color='IndianRed',linestyle='dashed', label='Shark v1.1 (L18)')
    ax.fill_between(xplot,yplot,yplot-errdn, facecolor='LightSalmon', linestyle='solid', alpha=0.4,interpolate=True)
    ax.fill_between(xplot,yplot,yplot+errup, facecolor='LightSalmon', linestyle='solid', alpha=0.4,interpolate=True)

    #Lange et al. (2016)
    m,r = common.load_observation(obsdir, 'SizesAndAM/rbulge_L16.dat', [0,1])
    m = np.log10(m)
    r = np.log10(r)
    ax.plot(m[0:50], r[0:50], linestyle='dotted',color='k')
    ax.plot(m[50:90], r[50:90], linestyle='dotted',color='k')
    ax.plot(m[90:len(r)], r[90:len(r)], linestyle='dotted',color='k', label='Lange+16 (bulges)')
    plot_profuse_sizes(ax, component=1, color='Crimson',include_ambig = True)
    plot_profuse_sizes(ax, component=1, color='Crimson',include_ambig = False)
    plot_profuse_sizes(ax, component=5, color='darkorange')

    common.prepare_legend(ax, ['r','IndianRed','black','Crimson','Crimson','darkorange'], loc=2)

    lm, lr, count, bs_err = common.load_observation(obsdir, 'SizeMass/GAMA_H-band_dlogM_0.25_reff.txt', [0,1,2,3])
    # Total ##################################
    xtit="$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit="$\\rm log_{10} (\\rm r_{\\star}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.2, 1.6

    ax = fig.add_subplot(313)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))

    ax.text(10.6, 0.,'total stellar',fontsize=12)

    #Predicted size-mass for disks
    ind = np.where(rcomb[0,0,:] != 0)
    xplot = xmf_obs[ind]
    yplot = rcomb[0,0,ind]
    errdn = rcomb[0,1,ind]
    errup = rcomb[0,2,ind]
    #np.save(os.path.join(outdir,'sizemass.npy'), np.array([xplot, yplot[0]]))
    ax.plot(xplot,yplot[0], color = 'ForestGreen',
            linewidth=2,label="Shark v2.0")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='ForestGreen', linestyle='solid', alpha=0.4,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='ForestGreen', linestyle='solid', alpha=0.4,interpolate=True)

    mL18, rL18, err16L18, err84L18 = common.load_observation(obsdir, 'Models/SharkVariations/SizeTotal_Lagos18.dat', [0,1,2,3])
    ind = np.where(rL18 != 0)
    xplot = mL18[ind]
    yplot = rL18[ind]
    errdn = err16L18[ind]
    errup = err84L18[ind]
    ax.plot(xplot,yplot,color='LimeGreen',linestyle='dashed', label='Shark v1.1 (L18)')
    ax.fill_between(xplot,errdn, errup, facecolor='LightGreen', linestyle='solid', alpha=0.4,interpolate=True)

    # Add GAMA H-band observations with bootstrapped error
    #ax.errorbar(lm, lr,yerr = [bs_err, bs_err], marker ='v',
    #         ls = 'none', mfc = 'None', markersize=5,
    #         color = 'gray', label = 'Lange+15 (H-band)')
    plot_profuse_sizes(ax, component=2, color='DarkSlateGray', include_ambig = True, r_bulge_assumed = 1)
    #plot_profuse_sizes(ax, component=2, color='DarkSlateGray', include_ambig = True, r_bulge_assumed = 0.1)


    common.prepare_legend(ax, ['ForestGreen','LimeGreen','DarkSlateGray'], loc=2)
    plt.tight_layout()
    common.savefig(outdir, fig, 'sizes_angular_momentum_model.pdf')


    #all disks and all bulges
    fig = plt.figure(figsize=(5,8.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star,disk}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm r_{\\star,disk}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.1, 1.5
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # LTG ##################################
    ax = fig.add_subplot(211)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(8.1,1.7,'disks all galaxies',fontsize=12)

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(disk_size[0,0,:] != 0)
    xplot = xmf_obs[ind]
    yplot = disk_size[0,0,ind]
    errdn = disk_size[0,1,ind]
    errup = disk_size[0,2,ind]
    ax.plot(xplot,yplot[0],color='b',linestyle='solid',label="Shark v2.0")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.3,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.3,interpolate=True)

    rdisk_am   = rb[30:59]
    rdisk_am16 = r16[30:59]
    rdisk_am84 = r84[30:59]
    ind = np.where(rdisk_am != 0)
    xplot = xmf[ind]
    yplot = rdisk_am[ind]
    errdn = rdisk_am16[ind]
    errup = rdisk_am84[ind]
    ax.plot(xplot,yplot,color='b',linestyle='dashed', label="Shark v1.1 (L18)")
    ax.fill_between(xplot,yplot,yplot-errdn, facecolor='b', linestyle='solid', alpha=0.8,interpolate=True)
    ax.fill_between(xplot,yplot,yplot+errup, facecolor='b', linestyle='solid', alpha=0.8,interpolate=True)

    #Lange et al. (2016)
    m,r = common.load_observation(obsdir, 'SizesAndAM/rdisk_L16.dat', [0,1])
    m = np.log10(m)
    r = np.log10(r)
    ax.plot(m[0:39], r[0:39], linestyle='dotted',color='k')
    ax.plot(m[39:78], r[39:78], linestyle='dotted',color='k')
    ax.plot(m[78:len(r)], r[78:len(r)], linestyle='dotted',color='k')

    common.prepare_legend(ax, ['b','b'], bbox_to_anchor=(0.005, 0.62))

    # ETGs ##################################
    xtit = "$\\rm log_{10} (\\rm M_{\\star,bulge}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm r_{\\star,bulge}/kpc)$"
    xmin, xmax, ymin, ymax = 8, 12, -0.35, 1.5
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    ax = fig.add_subplot(313)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(8.1,1.7,'bulges of all galaxies',fontsize=12)

    #Predicted size-mass for bulges in bulge-dominated systems
    ind = np.where(bulge_size[0,0,:] != 0)
    if(len(xmf[ind]) > 0):
        xplot = xmf_obs[ind]
        yplot = bulge_size[0,0,ind]
        errdn = bulge_size[0,1,ind]
        errup = bulge_size[0,2,ind]
        ax.plot(xplot,yplot[0],color='r',linestyle='solid')
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.3,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.3,interpolate=True)

    rdisk_am   = rb[60:99]
    rdisk_am16 = r16[60:99]
    rdisk_am84 = r84[60:99]

    ind = np.where(rdisk_am != 0)
    xplot = xmf[ind]
    yplot = rdisk_am[ind]
    errdn = rdisk_am16[ind]
    errup = rdisk_am84[ind]
    ax.plot(xplot,yplot,color='r',linestyle='dashed')
    ax.fill_between(xplot,yplot,yplot-errdn, facecolor='LightCoral', linestyle='solid', alpha=0.4,interpolate=True)
    ax.fill_between(xplot,yplot,yplot+errup, facecolor='LightCoral', linestyle='solid', alpha=0.4,interpolate=True)

    #Lange et al. (2016)
    m,r = common.load_observation(obsdir, 'SizesAndAM/rbulge_L16.dat', [0,1])
    m = np.log10(m)
    r = np.log10(r)
    ax.plot(m[0:50], r[0:50], linestyle='dotted',color='k')
    ax.plot(m[50:90], r[50:90], linestyle='dotted',color='k')
    ax.plot(m[90:len(r)], r[90:len(r)], linestyle='dotted',color='k')

    common.savefig(outdir, fig, 'sizes_angular_momentum_disks_bulges.pdf')



    fig = plt.figure(figsize=(9.5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm v_{\\rm max}/km s^{-1}$)"
    xmin, xmax, ymin, ymax = 10, 15, 1.3, 3.5
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (221, 222, 223, 224)
    indz = (0, 1, 2, 3)

    # choose type of selection:
    selec = 0 #all galaxies

    # LTG ##################################
    for z,s,p in zip(zlist, indz, subplots):
            ax = fig.add_subplot(p)
            common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
            ax.text(xleg, yleg, 'z=%s' % str(z))
            
            ind = np.where(vmax_halo_gal[s,0,:,selec] != 0)
            xplot = xmfh[ind]
            yplot = vmax_halo_gal[s,0,ind,selec]
            errdn = vmax_halo_gal[s,1,ind,selec]
            errup = vmax_halo_gal[s,2,ind,selec]
            ax.plot(xplot,yplot[0],color='k',label="central subhalos")
            ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.2,interpolate=True)
            ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.2,interpolate=True)

            common.prepare_legend(ax, ['k'], loc=2)


    common.savefig(outdir, fig, 'vmax_vs_subhalo.pdf')


def plot_specific_am(plt, outdir, obsdir, sam_stars_disk, sam_gas_disk_atom, sam_gas_disk_mol, sam_gas_disk_atom_ms, sam_halo, sam_bar, sam_stars, sam_stars2, sam_barv2, mass_cut, sam_gas_disk_atom2):

    fig = plt.figure(figsize=(4.5,4.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm j_{\\rm disk}/kpc\\, km s^{-1})$"
    xmin, xmax, ymin, ymax = 8, 11.5, 1.5, 5
    xleg = xmax - 0.5 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    # choose type of selection:
    selec = 1 #disk-dominated galaxies
    s = 0
    # LTG ##################################
    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'z=0 Shark-default', fontsize = 12)

    ind = np.where(sam_halo[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = sam_halo[s,0,ind,selec] + 3.0
    errdn = sam_halo[s,1,ind,selec]
    errup = sam_halo[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='k',label="DM")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.2,interpolate=True)

    #Predicted sAM-mass for disks in disk=dominated galaxies
    ind = np.where(sam_stars[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = sam_stars[s,0,ind,selec]+ 3.0
    errdn = sam_stars[s,1,ind,selec]
    errup = sam_stars[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='r',linestyle='dashed',label="stars")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.25,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.25,interpolate=True)

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(sam_gas_disk_atom[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = sam_gas_disk_atom[s,0,ind,selec] + 3.0
    errdn = sam_gas_disk_atom[s,1,ind,selec]
    errup = sam_gas_disk_atom[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='b',linestyle='dotted',label="atomic ISM")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.25,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.25,interpolate=True)

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(sam_gas_disk_mol[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = sam_gas_disk_mol[s,0,ind,selec] + 3.0
    errdn = sam_gas_disk_mol[s,1,ind,selec]
    errup = sam_gas_disk_mol[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='g',linestyle='dashdot',label="molecular ISM")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', alpha=0.25,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', alpha=0.25,interpolate=True)

    common.prepare_legend(ax, ['k'], loc=2)

    common.savefig(outdir, fig, 'specific_am_onlyz0.pdf')

    fig = plt.figure(figsize=(9.5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm j_{\\rm disk}/kpc\\, km s^{-1})$"
    xmin, xmax, ymin, ymax = 8, 11.5, 1.5, 5
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (221, 222, 223, 224)
    indz = (0, 1, 2, 3)

    # choose type of selection:
    selec = 1 #disk-dominated galaxies
    # LTG ##################################
    for z,s,p in zip(zlist, indz, subplots):
        ax = fig.add_subplot(p)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
        ax.text(xleg, yleg, 'z=%s' % str(z))
    
        ind = np.where(sam_halo[s,0,:,selec] != 0)
        xplot = xmf[ind]
        yplot = sam_halo[s,0,ind,selec] + 3.0
        errdn = sam_halo[s,1,ind,selec]
        errup = sam_halo[s,2,ind,selec]
        ax.plot(xplot,yplot[0],color='k',label="DM")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.2,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.2,interpolate=True)
    
        #Predicted sAM-mass for disks in disk=dominated galaxies
        ind = np.where(sam_stars[s,0,:,selec] != 0)
        xplot = xmf[ind]
        yplot = sam_stars[s,0,ind,selec]+ 3.0
        errdn = sam_stars[s,1,ind,selec]
        errup = sam_stars[s,2,ind,selec]
        ax.plot(xplot,yplot[0],color='r',label="stars")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.5,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.5,interpolate=True)
    
        #Predicted size-mass for disks in disk=dominated galaxies
        ind = np.where(sam_gas_disk_atom[s,0,:,selec] != 0)
        xplot = xmf[ind]
        yplot = sam_gas_disk_atom[s,0,ind,selec] + 3.0
        errdn = sam_gas_disk_atom[s,1,ind,selec]
        errup = sam_gas_disk_atom[s,2,ind,selec]
        ax.plot(xplot,yplot[0],color='b',label="atomic ISM")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.5,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.5,interpolate=True)
    
        #Predicted size-mass for disks in disk=dominated galaxies
        ind = np.where(sam_gas_disk_mol[s,0,:,selec] != 0)
        xplot = xmf[ind]
        yplot = sam_gas_disk_mol[s,0,ind,selec] + 3.0
        errdn = sam_gas_disk_mol[s,1,ind,selec]
        errup = sam_gas_disk_mol[s,2,ind,selec]
        ax.plot(xplot,yplot[0],color='g',label="molecular ISM")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', alpha=0.5,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', alpha=0.5,interpolate=True)
    
        common.prepare_legend(ax, ['k'], loc=2)


    common.savefig(outdir, fig, 'specific_am.pdf')
 

    def plot_observations_AM(ax, bar_type=0, color='r', mass_cut='8'):

        #Read observational data.
        mbO14, msO14, mHIO14, jbO14, jsO14, jgO14, jmolO14 = common.load_observation(obsdir, 'SizesAndAM/Obreschkow14_FP.dat', [6,7,8,11,12,13,15])

        ms_h22, mb_h22, js_h22, jb_h22   = common.load_observation(obsdir, 'SizesAndAM/Hardwick+22ab_SubsetForClaudia.dat', [4,8,10,12])

        ms, mserr, js, jserr  = common.load_observation(obsdir, 'SizesAndAM/Posti18_AMdata.dat', [3,5,6,8])
        ind = np.where(np.log10(ms) > mass_cut)
        if(len(ms[ind]) > 0):
           errmdn = np.log10(ms[ind]) - np.log10(ms[ind] - mserr[ind])
           errmup = np.log10(ms[ind] + mserr[ind]) - np.log10(ms[ind])
           errjdn = np.log10(js[ind]) - np.log10(js[ind] - jserr[ind])
           errjup = np.log10(js[ind] + jserr[ind]) - np.log10(js[ind]) 
           ms = ms[ind]
           js = js[ind] 
    
        mgB17, errmgB17, msB17, errmsB17, mbB17, errmbB17, jgB17, errjgB17, jsB17, errjsB17, jbB17, errjbB17 = common.load_observation(obsdir, 'SizesAndAM/LITTLETHINGS_Butler16.dat', [1,2,3,4,5,6,7,8,9,10,11,12])

        ind = np.where(msB17 > mass_cut)
        if(len(msB17[ind]) > 0):
           mgB17    = mgB17[ind]
           msB17    = msB17[ind]     
           errmsB17 = errmsB17[ind] 
           errmgB17 = errmgB17[ind]
           mbB17    = mbB17[ind]    
           errmbB17 = errmbB17[ind] 
           jgB17    = jgB17[ind]    
           errjgB17 = errjgB17[ind]  
           jsB17    = jsB17[ind]    
           errjsB17 = errjsB17[ind] 
           jbB17    = jbB17[ind]    
           errjbB17 = errjbB17[ind] 
    
        msD21, jsD21, errjsD21 = common.load_observation(obsdir, 'SizesAndAM/DiTeodoro21.dat', [0,3,4])
        ind = np.where(msD21 > mass_cut)
        if(len(msD21[ind]) > 0):
           msD21 = msD21[ind]
           jsD21 = jsD21[ind]
           errjsD21 = errjsD21[ind]
 
        msC17, mgC17, mbC17, errupmbC17, errdnmbC17, JstarC17, JgasC17, errupJgasC17, errdnJgasC17, jbarC17, errupjbarC17, errdnjbarC17 = common.load_observation(obsdir, 'SizesAndAM/Chowdhury17.dat', [1,2,5,6,7,8,8,10,11,15,16,17])
        ind = np.where(msC17 > mass_cut)
        if(len(msC17[ind]) > 0):
           msC17           = msC17[ind]      
           mgC17           = mgC17[ind]       
           mbC17           = mbC17[ind]       
           errupmbC17      = errupmbC17[ind]  
           errdnmbC17      = errdnmbC17[ind]  
           JstarC17        = JstarC17[ind]    
           JgasC17         = JgasC17[ind]     
           errupJgasC17    = errupJgasC17[ind]
           errdnJgasC17    = errdnJgasC17[ind]
           jbarC17         = jbarC17[ind]     
           errupjbarC17    = errupjbarC17[ind]
           errdnjbarC17    = errdnjbarC17[ind]
    
        #mbE17, jbE17 = common.load_observation(obsdir, 'SizesAndAM/Elson17.dat', [5,7]) 
        #mbE17 = np.log10(mbE17 * 1e8)
        #jbE17 = np.log10(jbE17)
       
        #IDGalaxy Mstar errup errdn Mgas errup errdn Mbar errup errdn jstar errup errdn jgas errup errdn  jbar  errup errdn
        #(msK18, errupmsK18, errdnmsK18, mgK18, errupmgK18, errdnmgK18, mbK18, errupmbK18, errdnmbK18, jsK18, errupjsK18, errdnjsK18, jgK18, 
        #errupjgK18, errdnjgK18, jbarK18, errupjbarK18, errdnjbarK18) = common.load_observation(obsdir, 'SizesAndAM/Kurapati18.dat', [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]) 
    
        mHIW, mHIWerr, msW, msWerr, mbW, mbWerr, jbW, jbWerr, jHIW, jsW = common.load_observation(obsdir, 'SizesAndAM/WHISP_Mass_AM_values_CLU.dat', [0,1,2,3,4,5,6,7,8,9])
        ind = np.where(msW > mass_cut)
        if(len(msW[ind]) > 0):
           mHIW    = mHIW[ind]   
           mHIWerr = mHIWerr[ind]
           msW     = msW[ind]    
           msWerr  = msWerr[ind] 
           mbW     = mbW[ind]    
           mbWerr  = mbWerr[ind] 
           jbW     = jbW[ind]    
           jbWerr  = jbWerr[ind] 
           jHIW    = jHIW[ind]   
           jsW     = jsW[ind]    

        msMP, msMPerr, mgMP, mgMPerr, mbMP, mbMPerr, jsMP, jsMPerr,  jgMP, jgMPerr, jbMP, jbMPerr = common.load_observation(obsdir, 'SizesAndAM/ManceraPina20_AM.dat', [0,1,2,3,4,5,6,7,8,9,10,11])
        ind = np.where(np.log10(msMP) > mass_cut)
        if(len(msMP[ind]) > 0):
           msMP    = msMP[ind]   
           msMPerr = msMPerr[ind]
           mgMP    = mgMP[ind]   
           mgMPerr = mgMPerr[ind]
           mbMP    = mbMP[ind]   
           mbMPerr = mbMPerr[ind]
           jsMP    = jsMP[ind]   
           jsMPerr = jsMPerr[ind]
           jgMP    = jgMP[ind]   
           jgMPerr = jgMPerr[ind]
           jbMP    = jbMP[ind]   
           jbMPerr = jbMPerr[ind]

        if bar_type == 0:
           ax.errorbar(np.log10(ms), np.log10(js), yerr=[errjdn,errjup], xerr=[errmdn,errmup],ls='None', mfc='None', ecolor = color, mec= color,marker='+',label="Posti+18")
           ax.errorbar(msC17, JstarC17-msC17, xerr=0.2, yerr=0.1, marker='p', color=color, fillstyle='none',label="Chowdhury+17", ls='None')
           ax.plot(ms_h22, np.log10(js_h22), marker='.', color=color, fillstyle='none', label="Hardwick+22a,b",  ls='None')

           ax.errorbar(msB17, jsB17, yerr=[errjsB17,errjsB17], xerr=[errmsB17,errmsB17], ls='None', mfc='None', ecolor = color, mec= color ,marker='s',label="Butler+17")
           ax.errorbar(msD21, jsD21, yerr=errjsD21, ls='None', mfc='None', ecolor = color, mec=color ,marker='D',label="Di Teodoro+21")
           #ax.errorbar(msK18, jsK18, yerr=[errdnjsK18,errupjsK18], xerr=[errdnmsK18,errupmsK18],ls='None', mfc='None', ecolor = color, mec= color,marker='*',label="Kurapati+18")
           ax.errorbar(msW, jsW, xerr=msWerr, yerr=0.1, ls='None', color=color, marker='*', label='Murugeshan+19')
           ind = np.where(jsMP != 0)
           ax.errorbar(np.log10(msMP[ind]), np.log10(jsMP[ind]), xerr = np.log10(msMP[ind] + msMPerr[ind]) - np.log10(msMP[ind]), yerr=np.log10(jsMP[ind] + jsMPerr[ind]) - np.log10(jsMP[ind]), ls='None', mfc='None', ecolor = color, mec= color,marker='o',label="Mancera-Pina+20")
           ax.plot(msO14, jsO14, marker='s', color='red', fillstyle='full',  ls='None', label='OG14')

        elif bar_type == 2:
           ax.errorbar(mgB17, jgB17, yerr=[errjgB17,errjgB17], xerr=errmgB17,ls='None', mfc='None', ecolor = color, mec= color,marker='s')
           ##ax.plot(msO14, jgO14, marker='o', color=color, fillstyle='none',  ls='None')
           ax.errorbar(mgC17, JgasC17-mgC17, xerr = 0.2, yerr=[errdnJgasC17, errupJgasC17], marker='p', color=color, fillstyle='none',  ls='None')
           #ax.errorbar(mgK18, jgK18, yerr=[errdnjgK18,errupjgK18], xerr=[errdnmgK18, errupmgK18],ls='None', mfc='None', ecolor = color, mec= color,marker='*')
           ax.errorbar(mHIW + np.log10(XH), jHIW, xerr=mHIWerr, yerr=0.1, ls='None', color=color, marker='*')
           ind = np.where(jgMP != 0)
           ax.errorbar(np.log10(mgMP[ind]), np.log10(jgMP[ind]), xerr = np.log10(mgMP[ind] + mgMPerr[ind]) - np.log10(mgMP[ind]), yerr=np.log10(jgMP[ind] + jgMPerr[ind]) - np.log10(jgMP[ind]), ls='None', mfc='None', ecolor = color, mec= color,marker='o')
           ax.plot(mHIO14, jgO14, marker='s', color='blue', fillstyle='full',  ls='None')

        elif bar_type == 3:
           ax.plot(mb_h22, np.log10(jb_h22), marker='.', color=color, fillstyle='none', ls='None')
           ax.errorbar(mbB17, jbB17, yerr=[errjbB17,errjbB17], xerr=[errmbB17,errmbB17],ls='None', mfc='None', ecolor = color, mec= color,marker='s')
           #ax.plot(mbO14, jbO14, marker='o', color=color, fillstyle='none', ls='None')
           ax.errorbar(mbC17, jbarC17, yerr=[errdnjbarC17,errupjbarC17], xerr=[errdnmbC17,errupmbC17], ls='None', mfc='None', ecolor = color, mec= color, marker='p')
           #ax.errorbar(mbK18, jbarK18, yerr=[errdnjbarK18,errupjbarK18], xerr=[errdnmbK18,errupmbK18],ls='None', mfc='None', ecolor = color, mec=color,marker='*')
           #ax.plot(mbE17, jbE17, marker='^', color=color, fillstyle='none', label='Elson17', ls='None')
           ax.errorbar(mbW, jbW, xerr=mbWerr, yerr=jbWerr, ls='None', color=color, marker='*')
           ind = np.where(jbMP != 0)
           ax.errorbar(np.log10(mbMP[ind]), np.log10(jbMP[ind]), xerr = np.log10(mbMP[ind] + mbMPerr[ind]) - np.log10(mbMP[ind]), yerr=np.log10(jbMP[ind] + jbMPerr[ind]) - np.log10(jbMP[ind]), ls='None', mfc='None', ecolor = color, mec= color,marker='o')
           ax.plot(mbO14, jbO14, marker='s', color='black', fillstyle='full',  ls='None')



    #plot angular momentum components separately. 
    fig = plt.figure(figsize=(15,5))
    s = 0 
    #plot stars
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm j/kpc\\, km s^{-1})$"
    xmin, xmax, ymin, ymax = mass_cut, 11.5, 1.5, 4.2
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymin + 0.1 * (ymax - ymin)

    ax = fig.add_subplot(141)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, "Stars")


    plot_observations_AM(ax, bar_type=0, color='grey', mass_cut=mass_cut)
    jst, jmole, jatomic, jbar = common.load_observation(obsdir, 'Models/SharkVariations/AngularMomentum.dat', [0,1,2,3])
    jsL18    = np.zeros(shape = (3, len(xmf)))
    jmolL18  = np.zeros(shape = (3, len(xmf)))
    jatomL18 = np.zeros(shape = (3, len(xmf)))
    jbarL18  = np.zeros(shape = (3, len(xmf)))
    i = 0
    p = 0
    for j in range(0,int(len(jst)/2)):
       jsL18[i,p]    = jst[j]
       jmolL18[i,p]  = jmole[j]
       jatomL18[i,p] = jatomic[j]
       jbarL18[i,p]  = jbar[j]
       p = p + 1
       if(p >= len(xmf)):
          p = 0
          i = i +1
   
    #Predicted sAM-mass for disks in disk=dominated galaxies
    ind = np.where(sam_stars[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = sam_stars[s,0,ind,selec]+ 3.0
    errdn = sam_stars[s,1,ind,selec]
    errup = sam_stars[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='r')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', alpha=0.35,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', alpha=0.35,interpolate=True)

    ind = np.where(sam_stars_disk[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = sam_stars_disk[s,0,ind,selec]+ 3.0
    errdn = sam_stars_disk[s,1,ind,selec]
    errup = sam_stars_disk[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='r', linestyle='dotted')

    ind = np.where(jsL18[0,:] != 0)
    xplot = xmf[ind]
    yplot = jsL18[0,ind]+ 3.0
    errdn = jsL18[1,ind]
    errup = jsL18[2,ind]
    ax.plot(xplot,yplot[0],color='r',linestyle='dashed')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='r', linestyle='dashed', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='r', linestyle='dashed', alpha=0.2,interpolate=True)

    common.prepare_legend(ax, ['grey', 'grey', 'grey','grey','grey','grey','grey'], loc=2)

    
    #plot molecular gas
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "" #$\\rm log_{10} (\\rm j_{\\rm H_2}/kpc\\, km s^{-1})$"
    ax = fig.add_subplot(142)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, "$\\rm H_{2}$")

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(sam_gas_disk_mol[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = sam_gas_disk_mol[s,0,ind,selec] + 3.0
    errdn = sam_gas_disk_mol[s,1,ind,selec]
    errup = sam_gas_disk_mol[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='g', label="ISM/stars AM transfer")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', alpha=0.35,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', alpha=0.35,interpolate=True)

    ind = np.where(jmolL18[0,:] != 0)
    xplot = xmf[ind]
    yplot = jmolL18[0,ind]+ 3.0
    errdn = jmolL18[1,ind]
    errup = jmolL18[2,ind]
    ax.plot(xplot,yplot[0],color='g',linestyle='dashed', label="Lagos+18")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', linestyle='dashed', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', linestyle='dashed', alpha=0.2,interpolate=True)

    #ax.plot(msO14, jmolO14, 'go',fillstyle='none')
    common.prepare_legend(ax, ['k'], loc=2)

    #plot atomic gas
    xtit = "$\\rm log_{10} (\\rm M_{\\rm atom}/M_{\odot})$"
    ytit = "" #$\\rm log_{10} (\\rm j_{\\rm HI}/kpc\\, km s^{-1})$"
    ax = fig.add_subplot(143)
    common.prepare_ax(ax, 8.0, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, "HI")

    plot_observations_AM(ax, bar_type=2, color='grey', mass_cut=mass_cut)
    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(sam_gas_disk_atom[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = sam_gas_disk_atom[s,0,ind,selec] + 3.0
    errdn = sam_gas_disk_atom[s,1,ind,selec]
    errup = sam_gas_disk_atom[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='b')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.35,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.35,interpolate=True)

    ind = np.where(jatomL18[0,:] != 0)
    xplot = xmf[ind]
    yplot = jatomL18[0,ind]+ 3.0
    errdn = jatomL18[1,ind]
    errup = jatomL18[2,ind]
    ax.plot(xplot,yplot[0],color='b',linestyle='dashed')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', linestyle='dashed', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', linestyle='dashed', alpha=0.2,interpolate=True)

    #plot total baryon
    xtit = "$\\rm log_{10} (\\rm M_{\\rm bar}/M_{\odot})$"
    ytit = "" #$\\rm log_{10} (\\rm j_{\\rm bar}/kpc\\, km s^{-1})$"
    ax = fig.add_subplot(144)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg-0.4, yleg, "All Baryons")
    plot_observations_AM(ax, bar_type=3, color='grey', mass_cut=mass_cut)

    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(sam_bar[s,0,:,selec] != 0)
    xplot = xmf[ind]
    yplot = sam_bar[s,0,ind,selec] + 3.0
    errdn = sam_bar[s,1,ind,selec]
    errup = sam_bar[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='k')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.25,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.25,interpolate=True)

    ind = np.where(jbarL18[0,:] != 0)
    xplot = xmf[ind]
    yplot = jbarL18[0,ind]+ 3.0
    errdn = jbarL18[1,ind]
    errup = jbarL18[2,ind]
    ax.plot(xplot,yplot[0],color='k',linestyle='dashed')
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', linestyle='dashed', alpha=0.1,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', linestyle='dashed', alpha=0.1,interpolate=True)

    #common.prepare_legend(ax, ['k'], loc=2)
    plt.tight_layout()
    common.savefig(outdir, fig, 'specific_am_z0_components.pdf')


    #plot angular momentum components separately for different B/T selections
    fig = plt.figure(figsize=(15,6))
    #plot stars
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm j/kpc\\, km s^{-1})$"
    xmin, xmax, ymin, ymax = 8.65, 11.5, 1.8, 4.4
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymin + 0.1 * (ymax - ymin)

    ax = fig.add_subplot(141)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg-0.45, yleg, "Stars (v2.0)", fontsize=12)

    linestyles = ['dotted', 'solid', 'dashed', 'dashdot']
    linewidths = [4,2,2,2]
    labels = ['$\\rm D/T > 0$', '$\\rm D/T>  0.25$', '$\\rm D/T>  0.5$', '$\\rm D/T>  0.75$']

    plot_observations_AM(ax, bar_type=0, color='grey', mass_cut=mass_cut)
   
    #Predicted sAM-mass for disks in disk=dominated galaxies
    for c in range(0,4):
        ind = np.where(sam_stars[s,0,:,c] != 0)
        xplot = xmf[ind]
        yplot = sam_stars[s,0,ind,c]+ 3.0
        errdn = sam_stars[s,1,ind,c]
        errup = sam_stars[s,2,ind,c]
        ax.plot(xplot,yplot[0],color='r', linestyle=linestyles[c], linewidth = linewidths[c])
        if(c == 1):
           ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='r', alpha=0.35,interpolate=True)

    common.prepare_legend(ax, ['grey', 'red', 'grey', 'grey','grey','grey','grey','grey','grey'], loc=2)

    
    #ax.plot(msO14, jmolO14, 'go',fillstyle='none')
    #common.prepare_legend(ax, ['k'], loc=2)

    #plot atomic gas
    xtit = "$\\rm log_{10} (\\rm M_{\\rm gas}/M_{\odot})$"
    ytit = "" #$\\rm log_{10} (\\rm j_{\\rm HI}/kpc\\, km s^{-1})$"
    ax = fig.add_subplot(142)
    common.prepare_ax(ax, 8, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    #ax.text(xleg-0.4, yleg, "Atomic Gas")

    plot_observations_AM(ax, bar_type=2, color='grey', mass_cut=mass_cut)
    #Predicted size-mass for disks in disk=dominated galaxies
    for c in range(0,4):
        ind = np.where((sam_gas_disk_atom[s,0,:,c] != 0) & (xmf > 8))
        xplot = xmf[ind]
        yplot = sam_gas_disk_atom[s,0,ind,c] + 3.0
        errdn = sam_gas_disk_atom[s,1,ind,c]
        errup = sam_gas_disk_atom[s,2,ind,c]
        ax.plot(xplot,yplot[0],color='b', linestyle=linestyles[c], linewidth = linewidths[c])
        if(c == 1):
            ax.plot(xplot,yplot[0],color='b', linestyle=linestyles[c], label='atomic gas (v2.0)')
            ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='b', alpha=0.35,interpolate=True)
        ind = np.where((sam_gas_disk_atom2[s,0,:,c] != 0) & (xmf > 8))
        xplot = xmf[ind]
        yplot = sam_gas_disk_atom2[s,0,ind,c] + 3.0
        ax.plot(xplot,yplot[0],color='LightSeaGreen', linestyle=linestyles[c])
        if(c == 1):
            ax.plot(xplot,yplot[0],color='LightSeaGreen', linestyle=linestyles[c], label='all ISM gas (v2.0)')
    common.prepare_legend(ax, ['b','LightSeaGreen','k','k'], loc=4)


    #plot total baryon
    xtit = "$\\rm log_{10} (\\rm M_{\\rm bar}/M_{\odot})$"
    ytit = "" #$\\rm log_{10} (\\rm j_{\\rm bar}/kpc\\, km s^{-1})$"
    ax = fig.add_subplot(143)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg-0.85, yleg, "All Baryons (v2.0)", fontsize=12)
    plot_observations_AM(ax, bar_type=3, color='grey', mass_cut=mass_cut)

    #Predicted size-mass for disks in disk=dominated galaxies
    for c in range(0,4):
        ind = np.where(sam_barv2[s,0,:,c] != 0)
        xplot = xmf[ind]
        yplot = sam_barv2[s,0,ind,c] + 3.0
        errdn = sam_barv2[s,1,ind,c]
        errup = sam_barv2[s,2,ind,c]
        ax.plot(xplot,yplot[0],color='k', linestyle=linestyles[c], label=labels[c], linewidth = linewidths[c])
        if(c == 1):
           ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='k', alpha=0.25,interpolate=True)

    common.prepare_legend(ax, ['k','k','k','k'], loc=2)

    #plot molecular gas
    xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot})$"
    ytit = "" #$\\rm log_{10} (\\rm j_{\\rm H_2}/kpc\\, km s^{-1})$"
    ax = fig.add_subplot(144)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(10.7, yleg, "$\\rm D/T > 0.25$", fontsize=12)

    c = 1
    #Predicted size-mass for disks in disk=dominated galaxies
    ind = np.where(sam_stars[s,0,:,c] != 0)
    xplot = xmf[ind]
    yplot = sam_stars[s,0,ind,c]+ 3.0
    errdn = sam_stars[s,1,ind,c]
    errup = sam_stars[s,2,ind,c]
    ax.plot(xplot,yplot[0],color='r', linestyle='solid', label = 'Stars (v2.0)')
    ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='r', alpha=0.35,interpolate=True)

    ind = np.where(sam_gas_disk_mol[s,0,:,c] != 0)
    xplot = xmf[ind]
    yplot = sam_gas_disk_mol[s,0,ind,c] + 3.0
    errdn = sam_gas_disk_mol[s,1,ind,c]
    errup = sam_gas_disk_mol[s,2,ind,c]
    ax.plot(xplot,yplot[0],color='darkgreen', linestyle='solid', label = 'Molecular gas (v2.0)')
    ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='g', alpha=0.35,interpolate=True)

    ind = np.where(sam_gas_disk_atom_ms[s,0,:,c] != 0)
    xplot = xmf[ind]
    yplot = sam_gas_disk_atom_ms[s,0,ind,c] + 3.0
    errdn = sam_gas_disk_atom_ms[s,1,ind,c]
    errup = sam_gas_disk_atom_ms[s,2,ind,c]
    ax.plot(xplot,yplot[0],color='b', linestyle='solid', label = 'Atomic gas (v2.0)')
    ax.fill_between(xplot,yplot[0]+errup[0],yplot[0]-errdn[0], facecolor='b', alpha=0.25,interpolate=True)


    msL18, jsL18, jgL18 = common.load_observation(obsdir, 'Models/SharkVariations/SpecificAM_Lagos18.dat', [0,1,2])
    ax.plot(msL18, jsL18, color='r', linestyle='dotted', linewidth = 3, label = 'Stars (L18)')
    ax.plot(msL18, jgL18, color='b', linestyle='dotted', linewidth = 3, label = 'Atomic gas (L18)')

    msO14, jsO14, jgO14, jmolO14 = common.load_observation(obsdir, 'SizesAndAM/Obreschkow14_FP.dat', [7,12,14,15])
    ax.plot(msO14, jsO14, marker='s', color='red', fillstyle='full', alpha=0.5, ls='None', label='OG14 $j_{\\star}$')
    ax.plot(msO14, jmolO14, marker='s', color='green', fillstyle='full', alpha=0.5, ls='None', label='OG14 $j_{\\rm H_2}$')
    ax.plot(msO14, jgO14, marker='s', color='b', fillstyle='full', alpha=0.5, ls='None', label = 'OG14 $j_{\\rm HI}$')
    print("median differences OG14", np.median(jsO14-jmolO14), np.median(jsO14-jgO14), np.median(jmolO14-jgO14))
    common.prepare_legend(ax, ['red','darkgreen','blue','red','blue','red','green','blue'], loc=2)

    plt.tight_layout()
    common.savefig(outdir, fig, 'specific_am_z0_components_BTselec.pdf')


    #for c in range (0,2):
    #   print 'will change selection'
    #   for i in range (0,3):
    #        print 'will change within the same sample'
    #        for x,y,z,a in zip(sam_stars_disk[s,i,:,c],sam_gas_disk_mol[s,i,:,c],sam_gas_disk_atom[s,i,:,c],sam_bar[s,i,:,c]):
    #             print x,y,z,a


def plot_specific_am_ratio(plt, outdir, obsdir, sam_ratio_halo_disk, sam_ratio_halo_gal, sam_ratio_halo_disk_gas, 
                           sam_vs_sam_halo_disk, sam_vs_sam_halo_gal, sam_vs_sam_halo_disk_gas):

    fig = plt.figure(figsize=(9.5,9.5))
    xtit = "$\\rm log_{10} (\\rm M_{\\rm halo}/M_{\odot})$"
    ytit = "$\\rm log_{10} (\\rm j_{\\star}/j_{\\rm halo}$)"
    xmin, xmax, ymin, ymax = 10, 15, -3, 1
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (221, 222, 223, 224)
    indz = (0, 1, 2, 3)

    # choose type of selection:
    selec = 1 #disk-dominated galaxies

    # LTG ##################################
    for z,s,p in zip(zlist, indz, subplots):
        ax = fig.add_subplot(p)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
        ax.text(xleg, yleg, 'z=%s' % str(z))
    
        ind = np.where(sam_ratio_halo_gal[s,0,:,selec] != 0)
        xplot = xmfh[ind]
        yplot = sam_ratio_halo_gal[s,0,ind,selec]
        errdn = sam_ratio_halo_gal[s,1,ind,selec]
        errup = sam_ratio_halo_gal[s,2,ind,selec]
        ax.plot(xplot,yplot[0],color='k',label="all stars")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.2,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.2,interpolate=True)
    
        ind = np.where(sam_ratio_halo_disk[s,0,:,selec] != 0)
        xplot = xmfh[ind]
        yplot = sam_ratio_halo_disk[s,0,ind,selec]
        errdn = sam_ratio_halo_disk[s,1,ind,selec]
        errup = sam_ratio_halo_disk[s,2,ind,selec]
        ax.plot(xplot,yplot[0],color='g',label="disk stars")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', alpha=0.2,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', alpha=0.2,interpolate=True)
    
        ind = np.where(sam_ratio_halo_disk_gas[s,0,:,selec] != 0)
        xplot = xmfh[ind]
        yplot = sam_ratio_halo_disk_gas[s,0,ind,selec]
        errdn = sam_ratio_halo_disk_gas[s,1,ind,selec]
        errup = sam_ratio_halo_disk_gas[s,2,ind,selec]
        ax.plot(xplot,yplot[0],color='b',label="disk gas")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.2,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.2,interpolate=True)
    
        common.prepare_legend(ax, ['k'], loc=2)


    common.savefig(outdir, fig, 'specific_am_ratio.pdf')

    selec = 0 #disk-dominated galaxies

    #plot specific AM vs. specific AM 
    fig = plt.figure(figsize=(9,8))
    xtit = "$\\rm log_{10} (\\rm j_{\\rm halo}/kpc\,km\,s^{-1})$"
    ytit = "$\\rm log_{10} (\\rm j_{\\star},j_{\\star,disk},j_{\\rm gas,disk}/kpc\,km\,s^{-1}$)"
    xmin, xmax, ymin, ymax = 1,6,1,6
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    subplots = (221, 222, 223, 224)
    indz = (0, 1, 2, 3)
    zinplot = (0, 0.5, 1, 2) 

    # LTG ##################################
    for z,s,p in zip(zinplot, indz, subplots):
        ax = fig.add_subplot(p)
        common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
        ax.text(xleg, yleg, 'z=%s' % str(z), fontsize=10)
        #if(s == 0):
        #    ax.text(5.5,6.4,'Lagos+18',fontsize=14)
    
        ind = np.where(sam_vs_sam_halo_gal[s,0,:,selec] != 0)
        xplot = xlf[ind]+3
        yplot = sam_vs_sam_halo_gal[s,0,ind,selec]+3
        errdn = sam_vs_sam_halo_gal[s,1,ind,selec]
        errup = sam_vs_sam_halo_gal[s,2,ind,selec]
        ax.plot(xplot,yplot[0],color='k',label="all stars")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.2,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.2,interpolate=True)
    
        ind = np.where(sam_vs_sam_halo_disk[s,0,:,selec] != 0)
        xplot = xlf[ind]+3
        yplot = sam_vs_sam_halo_disk[s,0,ind,selec]+3
        errdn = sam_vs_sam_halo_disk[s,1,ind,selec]
        errup = sam_vs_sam_halo_disk[s,2,ind,selec]
        ax.plot(xplot,yplot[0],color='g',label="disk stars")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', alpha=0.2,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', alpha=0.2,interpolate=True)
    
        ind = np.where(sam_vs_sam_halo_disk_gas[s,0,:,selec] != 0)
        xplot = xlf[ind]+3
        yplot = sam_vs_sam_halo_disk_gas[s,0,ind,selec]+3
        errdn = sam_vs_sam_halo_disk_gas[s,1,ind,selec]
        errup = sam_vs_sam_halo_disk_gas[s,2,ind,selec]
        ax.plot(xplot,yplot[0],color='b',label="disk gas")
        ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.2,interpolate=True)
        ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.2,interpolate=True)
    
        xplot = [1,5]
        ax.plot(xplot,xplot,color='grey',linestyle='dotted')
        if(s == 0):
           common.prepare_legend(ax, ['k','g','b'], loc=2)

    common.savefig(outdir, fig, 'specific_am_halo_vs_galaxy.pdf')
 
    selec = 0 #all galaxies

    #plot specific AM vs. specific AM 
    fig = plt.figure(figsize=(5,5))
    xtit = "$\\rm log_{10} (\\rm j_{\\rm halo}/kpc\,km\,s^{-1})$"
    ytit = "$\\rm log_{10} (\\rm j_{\\star},j_{\\star,disk},j_{\\rm gas,disk}/kpc\,km\,s^{-1}$)"
    xmin, xmax, ymin, ymax = 1,6,1,6
    xleg = xmax - 0.2 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    s = 0
    # LTG ##################################
    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'z=0', fontsize=10)
    #if(s == 0):
    #    ax.text(5.5,6.4,'Lagos+18',fontsize=14)
    
    ind = np.where(sam_vs_sam_halo_gal[s,0,:,selec] != 0)
    xplot = xlf[ind]+3
    yplot = sam_vs_sam_halo_gal[s,0,ind,selec]+3
    errdn = sam_vs_sam_halo_gal[s,1,ind,selec]
    errup = sam_vs_sam_halo_gal[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='k',label="all stars")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='k', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='k', alpha=0.2,interpolate=True)

    ind = np.where(sam_vs_sam_halo_disk[s,0,:,selec] != 0)
    xplot = xlf[ind]+3
    yplot = sam_vs_sam_halo_disk[s,0,ind,selec]+3
    errdn = sam_vs_sam_halo_disk[s,1,ind,selec]
    errup = sam_vs_sam_halo_disk[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='g',label="disk stars")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='g', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='g', alpha=0.2,interpolate=True)

    ind = np.where(sam_vs_sam_halo_disk_gas[s,0,:,selec] != 0)
    xplot = xlf[ind]+3
    yplot = sam_vs_sam_halo_disk_gas[s,0,ind,selec]+3
    errdn = sam_vs_sam_halo_disk_gas[s,1,ind,selec]
    errup = sam_vs_sam_halo_disk_gas[s,2,ind,selec]
    ax.plot(xplot,yplot[0],color='b',label="disk gas")
    ax.fill_between(xplot,yplot[0],yplot[0]-errdn[0], facecolor='b', alpha=0.2,interpolate=True)
    ax.fill_between(xplot,yplot[0],yplot[0]+errup[0], facecolor='b', alpha=0.2,interpolate=True)

    xplot = [1,5]
    ax.plot(xplot,xplot,color='grey',linestyle='dotted')
    common.prepare_legend(ax, ['k','g','b'], loc=2)

    common.savefig(outdir, fig, 'specific_am_halo_vs_galaxy_z0.pdf')
 

def plot_lambda(plt, outdir, obsdir, lambdaH,  lambda_jiang, lambda_mass, bt, ms, ssfr_z0):

    lambda_allstar = lambda_jiang[3,:] 
    lambda_disk= lambda_jiang[0,:]
    lambda_star= lambda_jiang[1,:]
    lambda_gas = lambda_jiang[2,:]

    fig = plt.figure(figsize=(5,12))
    xtit = ""
    ytit = "$\\rm log_{10}(\\lambda_{\\star})$"
    xmin, xmax, ymin, ymax = -3, 0, -4, 0
    xleg = xmin + 0.02 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)
    cutms = 7.5

    med = np.zeros(shape = (3, len(xlf)))

    ax = fig.add_subplot(411)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'all galaxies, stars+gas')

    ind = np.where((lambdaH > 0) & (lambda_allstar > 0) & (lambda_allstar < 10) & (ms > cutms))
    if(len(lambdaH[ind]) > 0):
       #xdata = np.log10(lambdaH[ind])
       #ydata = np.log10(lambda_allstar[ind])
       #us.density_contour(ax, xdata, ydata, 30, 30) #, **contour_kwargs)

       coeff = np.corrcoef(np.log10(lambdaH[ind]),np.log10(lambda_allstar[ind]))
       ax.text(xmin + 0.02 * (xmax - xmin), ymax - 0.25 * (ymax - ymin), 'R=%s' % str(np.around(coeff[0,1], decimals=3)))

       med[:] = us.wmedians(x=np.log10(lambdaH[ind]), y=np.log10(lambda_allstar[ind]), xbins = xlf, low_numbers=True)
       ind    = np.where((med[0,:] != 0) & (med[0,:] > -5))
       xobs   = xlf[ind]
       yobs   = med[0,ind]
       yerrdn = med[1,ind]
       yerrup = med[2,ind]
       ax.errorbar(xobs, yobs[0], yerr=[yerrdn[0],yerrup[0]], mfc='None', ecolor = 'grey', mec='grey',linestyle='solid', color='k')
 
    #ind = np.where((lambdaH > 0) & (lambda_allstar > 0) & (lambda_allstar < 10) & (ms > 10))
    #med[:] = us.wmedians(x=np.log10(lambdaH[ind]), y=np.log10(lambda_allstar[ind]), xbins = xlf, low_numbers=True)
    #ind    = np.where((med[0,:] != 0) & (med[0,:] > -5))
    #xobs   = xlf[ind]
    #yobs   = med[0,ind]
    #ax.plot(xobs, yobs[0], color='k',linestyle='dotted', label='$\\rm log_{10}(M_{\\star}/M_{\\odot}) > 10$')

    #ind = np.where((lambdaH > 0) & (lambda_allstar > 0) & (lambda_allstar < 10) & (ms > 7.5) & (ms < 9))
    #med[:] = us.wmedians(x=np.log10(lambdaH[ind]), y=np.log10(lambda_allstar[ind]), xbins = xlf, low_numbers=True)
    #ind    = np.where((med[0,:] != 0) & (med[0,:] > -5))
    #xobs   = xlf[ind]
    #yobs   = med[0,ind]
    #ax.plot(xobs, yobs[0], color='k',linestyle='dashed', label='$\\rm 7.5<log_{10}(M_{\\star}/M_{\\odot})<9$')

    ax = fig.add_subplot(412)
    ytit = "$\\rm log_{10}(\\lambda_{\\rm disk})$"
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'disk-dominated, stars+gas')

    ind = np.where((lambdaH > 0) & (lambda_disk > 0) & (bt < 0.5) & (ms > cutms))
    #xdata = np.log10(lambdaH[ind])
    #ydata = np.log10(lambda_disk[ind])
    #us.density_contour(ax, xdata, ydata, 30, 30) #, **contour_kwargs)
    coeff = np.corrcoef(np.log10(lambdaH[ind]),np.log10(lambda_disk[ind]))
    ax.text(xmin + 0.02 * (xmax - xmin), ymax - 0.25 * (ymax - ymin), 'R=%s' % str(np.around(coeff[0,1], decimals=3)))

    med[:] = us.wmedians(x=np.log10(lambdaH[ind]), y=np.log10(lambda_disk[ind]), xbins = xlf, low_numbers=True)
    ind    = np.where((med[0,:] != 0) & (med[0,:] > -5))
    xobs   = xlf[ind]
    yobs   = med[0,ind]
    yerrdn = med[1,ind]
    yerrup = med[2,ind]
    ax.errorbar(xobs, yobs[0], yerr=[yerrdn[0],yerrup[0]], mfc='None', ecolor = 'grey', mec='grey',linestyle='solid', color='k')

    ax = fig.add_subplot(413)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xleg, yleg, 'disk-dominated, stars')

    ind = np.where((lambdaH > 0) & (lambda_star > 0) & (bt < 0.5) & (ms > cutms))
    xdata = np.log10(lambdaH[ind])
    ydata = np.log10(lambda_star[ind])
    us.density_contour(ax, xdata, ydata, 30, 30) #, **contour_kwargs)
    coeff = np.corrcoef(np.log10(lambdaH[ind]),np.log10(lambda_star[ind]))
    ax.text(xmin + 0.02 * (xmax - xmin), ymax - 0.25 * (ymax - ymin), 'R=%s' % str(np.around(coeff[0,1], decimals=3)))

    med[:] = us.wmedians(x=np.log10(lambdaH[ind]), y=np.log10(lambda_star[ind]), xbins = xlf, low_numbers=True)
    ind    = np.where((med[0,:] != 0) & (med[0,:] > -5))
    xobs   = xlf[ind]
    yobs   = med[0,ind]
    yerrdn = med[1,ind]
    yerrup = med[2,ind]
    ax.errorbar(xobs, yobs[0], yerr=[yerrdn[0],yerrup[0]], mfc='None', ecolor = 'grey', mec='grey',linestyle='solid', color='k')

    ax = fig.add_subplot(414)
    xtit = "$\\rm log_{10}(\\lambda_{\\rm halo})$"
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
    ax.text(xmin + 0.02 * (xmax - xmin), yleg, 'disk-dominated, gas')

    ind = np.where((lambdaH > 0) & (lambda_gas > 0) & (bt < 0.5) & (ms > cutms))
    #xdata = np.log10(lambdaH[ind])
    #ydata = np.log10(lambda_gas[ind])
    #us.density_contour(ax, xdata, ydata, 30, 30) #, **contour_kwargs)
    coeff = np.corrcoef(np.log10(lambdaH[ind]),np.log10(lambda_gas[ind]))
    ax.text(xmin + 0.02 * (xmax - xmin), ymax - 0.25 * (ymax - ymin), 'R=%s' % str(np.around(coeff[0,1], decimals=3)))

    med[:] = us.wmedians(x=np.log10(lambdaH[ind]), y=np.log10(lambda_gas[ind]), xbins = xlf, low_numbers=True)
    ind    = np.where((med[0,:] != 0) & (med[0,:] > -5))
    xobs   = xlf[ind]
    yobs   = med[0,ind]
    yerrdn = med[1,ind]
    yerrup = med[2,ind]
    ax.errorbar(xobs, yobs[0], yerr=[yerrdn[0],yerrup[0]], mfc='None', ecolor = 'grey', mec='grey',linestyle='solid', color='k')

    common.savefig(outdir, fig, 'lambda_relation.pdf')


def main(modeldir, outdir, redshift_table, subvols, obsdir):

    #file_hdf5_sed = "Shark-SED-eagle-rr14.hdf5"

    plt = common.load_matplotlib()
    fields = {'galaxies': ('mstars_disk', 'mstars_bulge', 'mstars_burst_mergers', 'mstars_burst_diskinstabilities',
                           'mstars_bulge_mergers_assembly', 'mstars_bulge_diskins_assembly', 'm_bh', 'rstar_disk', 'rstar_bulge', 
                           'rgas_disk', 'rgas_bulge','type', 
                           'specific_angular_momentum_disk_star', 'specific_angular_momentum_bulge_star',
                           'specific_angular_momentum_disk_gas', 'specific_angular_momentum_bulge_gas',
                           'specific_angular_momentum_disk_gas_atom', 'specific_angular_momentum_disk_gas_mol',
                           'lambda_subhalo', 'mvir_subhalo', 'mvir_hosthalo', 'matom_disk', 'mmol_disk', 'mgas_disk',
                           'matom_bulge', 'mmol_bulge', 'mgas_bulge','sfr_disk', 'sfr_burst','vmax_subhalo', 'l_x', 'l_y', 'l_z')}
    fields_sed = {'SED/ab_dust': ('bulge_d','bulge_m','bulge_t','disk','total'),}
    fields_sed_nod = {'SED/ab_nodust': ('bulge_d','bulge_m','bulge_t','disk','total')}

    # Loop over redshift and subvolumes
    sam_stars_disk    = np.zeros(shape = (len(zlist), 3, len(xmf), 4))
    sam_gas_disk_atom = np.zeros(shape = (len(zlist), 3, len(xmf), 4))
    sam_gas_disk_atom2= np.zeros(shape = (len(zlist), 3, len(xmf), 4))
    sam_gas_disk_mol  = np.zeros(shape = (len(zlist), 3, len(xmf), 4))
    sam_gas_disk_atom_ms = np.zeros(shape = (len(zlist), 3, len(xmf), 4))
    sam_bar           = np.zeros(shape = (len(zlist), 3, len(xmf), 4))
    sam_barv2         = np.zeros(shape = (len(zlist), 3, len(xmf), 4))
    sam_halo          = np.zeros(shape = (len(zlist), 3, len(xmf), 4))
    sam_stars         = np.zeros(shape = (len(zlist), 3, len(xmf), 4))
    sam_stars2        = np.zeros(shape = (len(zlist), 3, len(xmf), 4))

    sam_ratio_halo_disk  = np.zeros(shape = (len(zlist), 3, len(xmfh), 4))
    sam_ratio_halo_gal   = np.zeros(shape = (len(zlist), 3, len(xmfh), 4))
    sam_ratio_halo_disk_gas = np.zeros(shape = (len(zlist), 3, len(xmfh), 4))
    vmax_halo_gal        = np.zeros(shape = (len(zlist), 3, len(xmfh), 4)) 

    sam_vs_sam_halo_disk  = np.zeros(shape = (len(zlist), 3, len(xlf), 4))
    sam_vs_sam_halo_gal   = np.zeros(shape = (len(zlist), 3, len(xlf), 4))
    sam_vs_sam_halo_disk_gas = np.zeros(shape = (len(zlist), 3, len(xlf), 4))

    disk_size_sat = np.zeros(shape = (len(zlist), 3, len(xmf)))
    disk_size_cen = np.zeros(shape = (len(zlist), 3, len(xmf))) 
    bulge_size    = np.zeros(shape = (len(zlist), 3, len(xmf_obs)))
    disk_size     = np.zeros(shape = (len(zlist), 3, len(xmf_obs))) 
    rcomb = np.zeros(shape = (len(zlist), 3, len(xmf_obs)))

    for index, snapshot in enumerate(redshift_table[zlist]):
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        #seds = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed, subvols, file_hdf5_sed)
        #seds_nod = common.read_photometry_data_variable_tau_screen(modeldir, snapshot, fields_sed_nod, subvols, file_hdf5_sed)
        #nbands = len(seds[0])

        (lh, lj, lm, bt, ms, ssfr, mass_cut)  = prepare_data(hdf5_data, index, sam_stars_disk, sam_gas_disk_atom, sam_gas_disk_atom2, sam_gas_disk_mol, sam_gas_disk_atom_ms, sam_halo, sam_ratio_halo_disk, 
                     sam_ratio_halo_gal, sam_ratio_halo_disk_gas, disk_size_sat, disk_size_cen, bulge_size, sam_vs_sam_halo_disk, sam_vs_sam_halo_gal,
                     sam_vs_sam_halo_disk_gas, sam_bar, sam_stars, sam_stars2,  vmax_halo_gal, sam_barv2, disk_size, rcomb) #, seds, seds_nod, nbands)

        if(index  == 0):
                lambdaH = lh
                lambda_jiang = lj
                lambda_mass  = lm
                BT_ratio = bt
                stellar_mass = ms
                ssfr_z0 = ssfr

    plot_specific_am(plt, outdir, obsdir, sam_stars_disk, sam_gas_disk_atom, sam_gas_disk_mol, sam_gas_disk_atom_ms, sam_halo, sam_bar, sam_stars, sam_stars2, sam_barv2, mass_cut, sam_gas_disk_atom2)
    #plot_specific_am_ratio(plt, outdir, obsdir, sam_ratio_halo_disk, sam_ratio_halo_gal, sam_ratio_halo_disk_gas, 
    #                       sam_vs_sam_halo_disk, sam_vs_sam_halo_gal, sam_vs_sam_halo_disk_gas)
    #plot_lambda(plt, outdir, obsdir, lambdaH, lambda_jiang, lambda_mass, BT_ratio, stellar_mass, ssfr_z0)
    plot_sizes(plt, outdir, obsdir, disk_size_cen, disk_size_sat, bulge_size, vmax_halo_gal, disk_size, rcomb)

if __name__ == '__main__':
    main(*common.parse_args())
