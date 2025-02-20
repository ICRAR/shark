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
"""QSO luminosity functions plot"""

import collections
import functools
import logging

import numpy as np
import matplotlib.cm as cm
import matplotlib.gridspec as gs
from matplotlib.lines import Line2D

import os
import common
import utilities_statistics as us

# Initialize arguments
zlist = (0.00, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10, 10.5, 11.0, 11.5, 12.0)
M_sun=1.9891e30
c_light=2.99792458e8
G=6.67259e-11
M_atom=1.66053873e-27
H_atom_mass=1.00794
sigma_Thomson=6.65245854e-29
J2erg=1e7
m2cm=1.0e2
kg2g=1.0e3
yr2s=3.1558e7
Lsun=3.839e33
eta_superEdd=4
alpha_ADAF=0.1
delta_ADAF=0.2
beta_ADAF=1-alpha_ADAF/0.55
ADAF_trans=1e-3*(delta_ADAF/5e-4)*((1-beta_ADAF)/beta_ADAF)*alpha_ADAF**2

print(f'The QSO luminosities are being calculated with alpha_ADAF={alpha_ADAF:.2f} and delta_ADAF={delta_ADAF:.2f}')

##################################
Llow=41.5
Lupp=47.5
Lbins=np.linspace(Llow,Lupp,16)
dL=Lbins[1]-Lbins[0]
Lx = Lbins + dL/2.0
max_logz=np.max(np.log10(1+np.array([zlist])))

mlow = 6
mupp = 12.5
dm = 0.125
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

slow = -3
supp = 3
dsfr = 0.125
sfrbins = np.arange(slow,supp,dsfr)
xsfr = sfrbins + dsfr/2.0


def load_lf_obs(obsdir):
    #Data from Hopkins+2007
    H07_LF=common.load_observation(obsdir,'AGN/Lbol_QLF_Hopkins2007.dat',[0,1,2,3,4])
    H07_LF[0,:]=np.where(H07_LF[0,:]==0.2,0.25,H07_LF[0,:])
    zsel=np.zeros(H07_LF.shape[1]).astype('bool')
    for z in zlist:
        zsel|=(H07_LF[0,:]==z).flatten()
    H07_LF=H07_LF[:,zsel]
    #Changing units to the same as Shark
    H07_LF[1,:]=3.9*10**(H07_LF[1,:]+33)
    H07_LF[2,:]=10**H07_LF[2,:]
    H07_LF[3,:]=10**H07_LF[3,:]
    
    #Data from Shankar+2009
    S09_LF=common.load_observation(obsdir,'AGN/Lbol_QLF_Shankar2009.dat',[0,1,2,3,4])
    S09_LF[0,:]=np.where(S09_LF[0,:]==0.1,0.0,S09_LF[0,:])
    S09_LF[0,:]=np.where(S09_LF[0,:]==1.1,1.0,S09_LF[0,:])
    S09_LF[0,:]=np.where(S09_LF[0,:]==2.1,2.0,S09_LF[0,:])
    S09_LF[0,:]=np.where(S09_LF[0,:]==3.8,4.0,S09_LF[0,:])
    S09_LF[2,:]=S09_LF[2,:]-S09_LF[1,:]
    S09_LF[3,:]=S09_LF[3,:]-S09_LF[1,:]
    S09_LF[4,:]=S09_LF[4,:]-S09_LF[1,:]
    S09_LF[1,:]=10**S09_LF[1,:]
    S09_LF[2,:]=10**S09_LF[2,:]
    S09_LF[3,:]=10**S09_LF[3,:]-S09_LF[2,:]
    S09_LF[4,:]=S09_LF[2,:]-10**S09_LF[4,:]
    
    #Data from Thorne+2022a
    T22_LF=common.load_observation(obsdir,'AGN/Lbol_QLF_Thorne2022a.dat',[0,1,2,3,4])
    T22_LF[3,:]=np.where(T22_LF[3,:]==0.02,0.00,T22_LF[3,:])
    T22_LF[3,:]=np.where(T22_LF[3,:]==0.20,0.25,T22_LF[3,:])
    T22_LF[3,:]=np.where(T22_LF[3,:]==0.45,0.50,T22_LF[3,:])
    T22_LF[3,:]=np.where(T22_LF[3,:]==1.75,2.00,T22_LF[3,:])
    zsel=np.zeros(T22_LF.shape[1]).astype('bool')
    for z in zlist:
        zsel|=T22_LF[3,:]==z
    T22_LF=T22_LF[:,zsel]
    T22_LF[0,:]=10**T22_LF[0,:]
    T22_LF[1,:]=10**T22_LF[1,:]
    T22_LF[2,:]=10**T22_LF[2,:]
    
    return(H07_LF,S09_LF,T22_LF)

def prepare_data(hdf5_data, snapshot, read_spin, index, densities_agn, sfrf, msf, LF_bol):

    if(read_spin):
       (h0, volh, MBH, bh_accretion_rate_hh, bh_accretion_rate_sb, BH_spin, sfrd, sfrb, mstarsd, mstarsb) = hdf5_data
    else:
        (h0, volh, MBH, bh_accretion_rate_hh, bh_accretion_rate_sb, sfrd, sfrb, mstarsd, mstarsb) = hdf5_data
        BH_spin = np.zeros(shape = len(MBH))
        BH_spin[:] = 0.68

    vol = volh/pow(h0,3.)

    MBH = MBH.astype('float64') / h0
    MBH_acc = (bh_accretion_rate_hh + bh_accretion_rate_sb) #.astype('float64')
    MBH_acc /= (h0 * 1e9)
    densities_agn[index,1] = sum(MBH_acc)/vol
    densities_agn[index,2] = sum(mstarsd + mstarsb)/h0/vol
    densities_agn[index,3] = sum(sfrd + sfrb)/h0/1e9/vol

    ind = np.where(mstarsd + mstarsb > 0)
    msf[index,:] = np.histogram(np.log10((mstarsd[ind] + mstarsb[ind])/h0),bins=np.append(mbins,mupp))[0]/vol/dm 
    ind = np.where(sfrd + sfrb > 0)
    sfrf[index,:] = np.histogram(np.log10((sfrd[ind] + sfrb[ind])/h0/1e9),bins=np.append(sfrbins, supp))[0]/vol/dsfr


    
    BH_spin = BH_spin.astype('float64')

    def acc_eff_calc(a):
        a1=np.abs(a)
        a2=a**2
        Z1=1+((1+a1)**(1/3.0)+(1-a1)**(1/3.0))*(1-a2)**(1/3.0)
        Z2=(3*a2+Z1**2)**0.5
        
        r_lso=3+Z2
        r_temp=((3-Z1)*(3+Z1+2*Z2))**0.5
        a_positive=a>=0
        r_lso[a_positive]-=r_temp[a_positive]
        r_lso[~a_positive]+=r_temp[~a_positive]
        
        acc_eff=1-(1-(2/(3*r_lso)))**0.5
        acc_eff=np.where(acc_eff<0,0.07,acc_eff)
        acc_eff=np.where(acc_eff>0.5,0.5,acc_eff)
        return(r_lso,acc_eff)
    
    #Efficiency
    r_lso,acc_eff=acc_eff_calc(BH_spin)

    #Eddington luminosity
    L_Edd=4*np.pi*c_light*G*M_sun*M_atom*H_atom_mass/(sigma_Thomson*1e40)
    L_Edd*=MBH*J2erg
   
    #Eddington MBH accretion rate and ratio
    M_dot_Edd=1e40*L_Edd/(0.1*(c_light*m2cm)**2)
    M_dot=MBH_acc*M_sun*kg2g/yr2s

    m_dot=np.where((M_dot_Edd>0)&(M_dot>0),M_dot/M_dot_Edd,np.nan)
    #Bolometric luminosity
    L_bol=np.zeros(len(m_dot))
    ADAF_low=(0<m_dot)&(m_dot<=ADAF_trans)
    ADAF_high=(m_dot>ADAF_trans)&(m_dot<1e-2)
    TD=m_dot>=1e-2
    
    temp=2e-4*acc_eff[ADAF_low]*M_dot[ADAF_low]*(c_light*m2cm)**2/(r_lso[ADAF_low]*1e40)
    L_bol[ADAF_low]=temp*6*(delta_ADAF/5e-4)*(1-beta_ADAF)/0.5
    
    temp=0.2*acc_eff[ADAF_high]*M_dot[ADAF_high]*(c_light*m2cm)**2/(r_lso[ADAF_high]*1e40)
    L_bol[ADAF_high]=temp*6*beta_ADAF*m_dot[ADAF_high]/(0.5*alpha_ADAF**2)
    
    L_bol[TD]=acc_eff[TD]*M_dot[TD]*(c_light*m2cm)**2/1e40
    
    #Correcting for super Eddington
    SE=TD&(m_dot>eta_superEdd*(0.1/acc_eff))
    L_bol[SE]=eta_superEdd*(1+np.log((m_dot[SE]/eta_superEdd)*(acc_eff[SE]/0.1)))*L_Edd[SE]
 
    L_bol=np.where(L_bol>0,np.log10(L_bol)+40,np.nan)

    ind=np.where(L_bol > 0)
    LF_bol[index,:] = np.histogram(L_bol[ind],bins=np.append(Lbins,Lupp))[0]/vol/dL
    densities_agn[index,0] = sum(10**L_bol[ind]/1e40)/vol
    print("Total bolometric density", densities_agn[index,0], densities_agn[index,1], " at redshift:", zlist[index])

def main(modeldir, outdir, redshift_table, subvols, obsdir):
    plt = common.load_matplotlib()
    
    densities_sm_agn = np.zeros(shape = (len(zlist), 4))
    LF_bol = np.zeros(shape = (len(zlist), len(Lbins)))
    msf = np.zeros(shape = (len(zlist), len(xmf)))
    sfrf = np.zeros(shape = (len(zlist), len(xsfr)))

    read_spin = True

    if(read_spin):
        fields = {'galaxies': ('m_bh', 'bh_accretion_rate_hh', 'bh_accretion_rate_sb', 'bh_spin', 'sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge')}
    else:
        fields = {'galaxies': ('m_bh', 'bh_accretion_rate_hh', 'bh_accretion_rate_sb', 'sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge')}
    
    for index, snapshot in enumerate(redshift_table[zlist]):
        print("will read data at z", zlist[index])
        hdf5_data = common.read_data(modeldir, snapshot, fields, subvols)
        prepare_data(hdf5_data, snapshot, read_spin, index, densities_sm_agn, sfrf, msf, LF_bol)
    
    for a,b,c,d,e in zip(zlist, densities_sm_agn[:,0], densities_sm_agn[:,1], densities_sm_agn[:,2], densities_sm_agn[:,3]):
        print(a,b,c,d,e)

    for index in range(0,len(zlist)):
        if(zlist[index] >= 5.5):
           print("#AGN bolometric LF at z:", zlist[index])
           for a,b in zip(Lx, LF_bol[index,:]):
               print(a,b)
    for index in range(0,len(zlist)):
        if(zlist[index] >= 5.5):
           print("#SMF at z:", zlist[index])
           for a,b in zip(xmf, msf[index,:]):
               print(a,b)

    for index in range(0,len(zlist)):
        if(zlist[index] >= 5.5):
           print("#SFR function z:", zlist[index])
           for a,b in zip(xsfr, sfrf[index,:]):
               print(a,b)

if __name__ == '__main__':
    main(*common.parse_args())

