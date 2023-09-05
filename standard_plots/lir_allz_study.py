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

zlist=np.array([0.194739, 0.450678, 0.8, 0.9, 1.20911, 1.59696, 2.00392, 2.47464723643932, 3.01916, 3.50099697082904, 3.95972, 4.465197621546, 5.02220991014863, 5.52950356184419, 5.96593,  6.55269895697227, 7.0, 8.0, 9.0])
#zlist=np.array([0.194739, 0.254144, 0.359789, 0.450678, 0.8, 0.849027, 0.9, 1.20911, 1.28174, 1.39519, 1.59696, 2.00392, 2.47464723643932, 2.76734390952347, 3.01916, 3.21899984389701, 3.50099697082904, 3.7248038025221, 3.95972, 4.465197621546, 4.73693842543988, 5.02220991014863, 5.2202206934302, 5.52950356184419, 5.74417977285603, 5.96593, 6.19496927748119, 6.55269895697227, 8.0, 9.0, 10.0])

#zlist=np.array([0.194739, 0.254144, 0.359789, 0.450678, 0.8, 0.849027, 0.9, 1.20911, 1.28174, 1.39519, 1.59696, 2.00392, 2.47464723643932, 2.76734390952347, 3.01916, 3.21899984389701, 3.50099697082904, 3.7248038025221, 3.95972, 4.465197621546, 4.73693842543988, 5.02220991014863, 5.2202206934302, 5.52950356184419, 5.74417977285603, 5.96593, 6.19496927748119, 6.55269895697227, 7.05756323172746, 7.45816170313544, 7.73629493731708, 8.02352,8.32018565809831, 8.47220854014322, 8.78358705435761, 8.94312532315157, 9.27010372804765, 9.437541750167, 9.78074128377067, 9.95655])

##################################
#Constants
RExp     = 1.67
MpcToKpc = 1e3
G        = 4.299e-9 #Gravity constant in units of (km/s)^2 * Mpc/Msun
c_light  = 299792458.0 #m/s
PI       = 3.141592654

mlow = 6.5
mupp = 14.0
dm = 0.2
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0


def prepare_data(hdf5_data, seds, index, LFs_dust, obsdir, nbrightgals, nbrightgals_err, csfr):

    (h0, volh, mdisk, mbulge, sfr_disk, sfr_burst) = hdf5_data

    sm_thresh = 3.98e10
    ssfr_thresh = 1e-10
    ssfr_thresh2 = 1e-11

    lir_total = seds[1] #total absolute magnitudes with dust
    mstar_tot = (mdisk + mbulge)/h0
    sfr_tot = (sfr_disk + sfr_burst)/h0/1e9
    csfr[index,1] = np.sum(sfr_tot)
    
    ind = np.where(mstar_tot > 0)
    mstar_tot = mstar_tot[ind]
    sfr_tot = sfr_tot[ind]
    lir_total = lir_total[0,:]
    csfr[index,0] = np.sum(lir_total) * 1.72e-10

    csfr[index,:] = np.log10(csfr[index,:] / (volh / h0**3))
    #no error added
    ind = np.where((lir_total >= 1e11) & (lir_total < 1e20))

    brightgals = lir_total[ind]
    nbrightgals[0,0,index] = brightgals.size
    nbrightgals[0,1,index] = sum(mstar_tot[ind])

    ind = np.where((lir_total >= 1e12) & (lir_total < 1e20))
    brightgals = lir_total[ind]
    nbrightgals[1,0,index] = brightgals.size
    nbrightgals[1,1,index] = sum(mstar_tot[ind])

    ind = np.where((lir_total >= 1e11) & (lir_total < 1e20) & (mstar_tot >= 1e11))
    brightgals = lir_total[ind]
    nbrightgals[2,0,index] = brightgals.size
    nbrightgals[2,1,index] = sum(mstar_tot[ind])

    ind = np.where(mstar_tot >= 1e11)
    brightgals = mstar_tot[ind]
    nbrightgals[3,0,index] = brightgals.size
    nbrightgals[3,1,index] = sum(mstar_tot[ind])

    ind = np.where((lir_total >= 1e11) & (lir_total < 1e20) & (mstar_tot >= sm_thresh))
    brightgals = lir_total[ind]
    nbrightgals[4,0,index] = brightgals.size
    nbrightgals[4,1,index] = sum(mstar_tot[ind])

    ind = np.where((lir_total >= 1e12) & (lir_total < 1e20) & (mstar_tot >= sm_thresh))
    brightgals = lir_total[ind]
    nbrightgals[8,0,index] = brightgals.size
    nbrightgals[8,1,index] = sum(mstar_tot[ind])

    ind = np.where(mstar_tot >= sm_thresh)
    brightgals = mstar_tot[ind]
    nbrightgals[5,0,index] = brightgals.size
    nbrightgals[5,1,index] = sum(mstar_tot[ind])

    ind = np.where((mstar_tot >= sm_thresh) & (sfr_tot/mstar_tot < ssfr_thresh)) 
    brightgals = mstar_tot[ind]
    nbrightgals[6,0,index] = brightgals.size
    nbrightgals[6,1,index] = sum(mstar_tot[ind])

    ind = np.where((mstar_tot >= 1e11) & (sfr_tot/mstar_tot < ssfr_thresh)) 
    brightgals = mstar_tot[ind]
    nbrightgals[7,0,index] = brightgals.size
    nbrightgals[7,1,index] = sum(mstar_tot[ind])

    ind = np.where((mstar_tot >= sm_thresh) & (sfr_tot/mstar_tot < ssfr_thresh2)) 
    brightgals = mstar_tot[ind]
    nbrightgals[9,0,index] = brightgals.size
    nbrightgals[9,1,index] = sum(mstar_tot[ind])

    ind = np.where((mstar_tot >= 1e11) & (sfr_tot/mstar_tot < ssfr_thresh2)) 
    brightgals = mstar_tot[ind]
    nbrightgals[10,0,index] = brightgals.size
    nbrightgals[10,1,index] = sum(mstar_tot[ind])


    mstar_tot = 10.0**(np.log10(mstar_tot) + np.random.normal(0,0.2,mstar_tot.size))

    ind = np.where(lir_total > 0)
    allgals = lir_total[ind]
    lir_total[ind] = 10.0**(np.log10(lir_total[ind]) + np.random.normal(0,0.2,allgals.size))

    #error added
    ind = np.where((lir_total >= 1e11) & (lir_total < 1e20))
    brightgals = lir_total[ind]
    nbrightgals_err[0,0,index] = brightgals.size
    nbrightgals_err[0,1,index] = sum(mstar_tot[ind])

    ind = np.where((lir_total >= 1e12) & (lir_total < 1e20))
    brightgals = lir_total[ind]
    nbrightgals_err[1,0,index] = brightgals.size
    nbrightgals_err[1,1,index] = sum(mstar_tot[ind])

    ind = np.where((lir_total >= 1e11) & (lir_total < 1e20) & (mstar_tot >= 1e11))
    brightgals = lir_total[ind]
    nbrightgals_err[2,0,index] = brightgals.size
    nbrightgals_err[2,1,index] = sum(mstar_tot[ind])

    ind = np.where(mstar_tot >= 1e11)
    brightgals = mstar_tot[ind]
    nbrightgals_err[3,0,index] = brightgals.size
    nbrightgals_err[3,1,index] = sum(mstar_tot[ind])

    ind = np.where((lir_total >= 1e11) & (lir_total < 1e20) & (mstar_tot >= sm_thresh))
    brightgals = lir_total[ind]
    nbrightgals_err[4,0,index] = brightgals.size
    nbrightgals_err[4,1,index] = sum(mstar_tot[ind])

    ind = np.where((lir_total >= 1e12) & (lir_total < 1e20) & (mstar_tot >= sm_thresh))
    brightgals = lir_total[ind]
    nbrightgals_err[8,0,index] = brightgals.size
    nbrightgals_err[8,1,index] = sum(mstar_tot[ind])

    ind = np.where(mstar_tot >= sm_thresh)
    brightgals = mstar_tot[ind]
    nbrightgals_err[5,0,index] = brightgals.size
    nbrightgals_err[5,1,index] = sum(mstar_tot[ind])

    ind = np.where((mstar_tot >= sm_thresh) & (sfr_tot/mstar_tot < ssfr_thresh)) 
    brightgals = mstar_tot[ind]
    nbrightgals_err[6,0,index] = brightgals.size
    nbrightgals_err[6,1,index] = sum(mstar_tot[ind])

    ind = np.where((mstar_tot >= 1e11) & (sfr_tot/mstar_tot < ssfr_thresh)) 
    brightgals = mstar_tot[ind]
    nbrightgals_err[7,0,index] = brightgals.size
    nbrightgals_err[7,1,index] = sum(mstar_tot[ind])

    ind = np.where((mstar_tot >= sm_thresh) & (sfr_tot/mstar_tot < ssfr_thresh2)) 
    brightgals = mstar_tot[ind]
    nbrightgals_err[9,0,index] = brightgals.size
    nbrightgals_err[9,1,index] = sum(mstar_tot[ind])

    ind = np.where((mstar_tot >= 1e11) & (sfr_tot/mstar_tot < ssfr_thresh2)) 
    brightgals = mstar_tot[ind]
    nbrightgals_err[10,0,index] = brightgals.size
    nbrightgals_err[10,1,index] = sum(mstar_tot[ind])



    #compute luminosity function
    ind = np.where((lir_total > 0) & (lir_total < 1e20))
    H, bins_edges = np.histogram(np.log10(lir_total[ind]),bins=np.append(mbins,mupp))
    LFs_dust[index,:] = LFs_dust[index,:] + H
    LFs_dust[index,:] = LFs_dust[index,:] / (volh/h0**3.0) / dm
    ind = np.where(LFs_dust[index,:] != 0)
    LFs_dust[index,ind] = np.log10(LFs_dust[index,ind])

    return(volh, h0)
    
def plot_lir_lf(plt, outdir, obsdir, LFs_dust, file_name, csfr):

    fig = plt.figure(figsize=(5,4))
    plt.subplots_adjust(left=0.15, bottom=0.17)
    xtit="$\\rm redshift$"
    ytit="$\\rm log_{10}(CSFR/M_{\\odot} yr^{-1} Mpc^{-3}]$"

    xmin, xmax, ymin, ymax = 0, 10.1, -4, -0.5

    ax = fig.add_subplot(111)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1))
    ax.plot(zlist, csfr[:,0], color='red', linestyle='solid', label='Shark LIR-SFR')
    ax.plot(zlist, csfr[:,0]-0.3, color='red', linestyle='dotted', label='Shark LIR-SFR -0.4dex')
    ax.plot(zlist, csfr[:,1], color='k', linestyle='solid', label='Shark total SFR')
    #for a,b,c in zip(zlist, csfr[:,0],csfr[:,1]):
    #    print(a,b,c)
    common.prepare_legend(ax, ['red','red','k'], loc='lower left')
    common.savefig(outdir, fig, "CSFR_"+file_name+".pdf")

    fig = plt.figure(figsize=(6,6))
    ytit = "$\\rm log_{10} (\\rm \\phi/\\, cMpc^{-3}\\, dex^{-1})$"
    xtit = "$\\rm log_{10} (L_{\\rm TIR}/L_{\\odot})$"
    xmin, xmax, ymin, ymax = 8.5, 12.8, -6, -1.5
    xleg = xmax - 0.3 * (xmax - xmin)
    yleg = ymax - 0.1 * (ymax - ymin)

    cols = ('Indigo','purple','Navy','MediumBlue','Green','MediumAquamarine','LightGreen','YellowGreen','Gold','Orange','Coral','OrangeRed','red','DarkRed','FireBrick','Crimson','IndianRed','LightCoral','Maroon','brown','Sienna','SaddleBrown','Chocolate','Peru','DarkGoldenrod','Goldenrod','SandyBrown')
    ax = fig.add_subplot(111)

    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.5, 0.5, 0.5, 0.5))

    for b in range(0,len(zlist)):
        inp = np.where(LFs_dust[b,:] != 0)
        x = xmf[inp]
        y = LFs_dust[b,inp]
        ax.plot(x, y[0], linestyle='solid',color=cols[b], label='z=%s' % str(zlist[b]))
    common.prepare_legend(ax, cols, loc=3)
    common.savefig(outdir, fig, 'LIR_total_allz_'+file_name+'.pdf')


def main(model_dir, outdir, redshift_table, subvols, obsdir):

    plt = common.load_matplotlib()

    file_name = "eagle-rr14"
    file_hdf5_sed = "Shark-SED-" + file_name + ".hdf5"
    fields_sed = {'SED/lir_dust': ('disk','total'),}

    fields = {'galaxies': ('mstars_disk', 'mstars_bulge','sfr_disk', 'sfr_burst')}

    #Bands information:
    #(0): "FUV_GALEX", "NUV_GALEX", "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS",
    #(6): "z_SDSS", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1_WISE",
    #(12): "I1_Spitzer", "I2_Spitzer", "W2_WISE", "I3_Spitzer", "I4_Spitzer",
    #(17): "W3_WISE", "W4_WISE", "P70_Herschel", "P100_Herschel",
    #(21): "P160_Herschel", "S250_Herschel", "S350_Herschel", "S450_JCMT",
    #(25): "S500_Herschel", "S850_JCMT", "Band9_ALMA", "Band8_ALMA",
    #(29): "Band7_ALMA", "Band6_ALMA", "Band5_ALMA", "Band4_ALMA"

    csfr = np.zeros(shape = (len(zlist), 2))
    LFs_dust     = np.zeros(shape = (len(zlist), len(mbins)))
    nbrightgals = np.zeros(shape = (11,2,len(zlist)))
    nbrightgals_err = np.zeros(shape = (11,2,len(zlist)))

    for index, snapshot in enumerate(redshift_table[zlist]):
        #print("Will read snapshot %s" % (str(snapshot)))
        hdf5_data = common.read_data(model_dir, snapshot, fields, subvols)
        seds = common.read_photometry_data_variable_tau_screen(model_dir, snapshot, fields_sed, subvols, file_hdf5_sed)
        (volh, h0) = prepare_data(hdf5_data, seds, index, LFs_dust, obsdir, nbrightgals, nbrightgals_err, csfr)

    nbrightgals = nbrightgals/(volh/h0**3.0)
    nbrightgals_err = nbrightgals_err/(volh/h0**3.0)

    ngals = nbrightgals
    j = 0


    Write_Tables = True 
    #for a,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19 in zip(xmf, LFs_dust[0,:], LFs_dust[1,:], LFs_dust[2,:], LFs_dust[3,:], LFs_dust[4,:], LFs_dust[5,:], LFs_dust[6,:], LFs_dust[7,:], LFs_dust[8,:], LFs_dust[9,:], LFs_dust[10,:], LFs_dust[11,:], LFs_dust[12,:], LFs_dust[13,:], LFs_dust[14,:], LFs_dust[15,:], LFs_dust[16,:], LFs_dust[17,:], LFs_dust[18,:]):
    #        print(a,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19)

    if(Write_Tables == True): 

       print('#Computed using the Shark model (Lagos et al. 2018; MNRAS.481.3573) using the SED modelling introduced in Lagos et al. (2019; MNRAS.489.4196). Also refer to Lagos et al. (2020; MNRAS.499.1948) for a detailed analysis of dusty galaxies in Shark.')
       print('#Table details:')
       print('#Number densities in units of Mpc^-3')
       print('#N1: number density of galaxies with LIR>1e11Lsun')
       print('#N2: number density of galaxies with LIR>1e12Lsun')
       print('#N3: number density of galaxies with LIR>1e11Lsun and Mstar>1e11Msun')
       print('#N4: number density of galaxies with Mstar>1e11Msun:')
       print('#N5: number density of galaxies with LIR>1e11Lsun and Mstar>3.9e10Msun')
       print('#N6: number density of galaxies with Mstar>3.9e10Msun')
       print('#N7: number density of galaxies with Mstar>3.9e10Msun & SSFR<1e-10yr^-1')
       print('#N8: number density of galaxies with Mstar>1e11Msun & SSFR<1e-10yr^-1')
       print('#N9: number density of galaxies with LIR>1e12Lsun and Mstar>3.9e10Msun')
       print('#N10: number density of galaxies with Mstar>3.9e10Msun & SSFR<1e-11yr^-1')
       print('#N11: number density of galaxies with Mstar>1e11Msun & SSFR<1e-11yr^-1')
       print('#(WE): refers to the cuts above applied after convolving IR luminosities and stellar masses with a gaussian of width 0.2dex')
       print('#SMD: refers to stellar mass density in units of Msun/Mpc^-3 and the sample numbers are as for the number densities.')
       print('# ')
       print('#columns:')
       print('#redshift N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11')
       for a,b,c,d,e,f,g,h,i,k,l,m in zip(zlist,ngals[0,j,:],ngals[1,j,:],ngals[2,j,:],ngals[3,j,:],ngals[4,j,:],ngals[5,j,:],ngals[6,j,:],ngals[7,j,:],ngals[8,j,:],ngals[9,j,:],ngals[10,j,:]):
           print(a,b,c,d,e,f,g,h,i,k,l,m)
       print('# ')
       print('#redshift N1(WE) N2(WE) N3(WE) N4(WE) N5(WE) N6(WE) N7(WE) N8(WE) N9(WE)')
       ngals = nbrightgals_err
       for a,b,c,d,e,f,g,h,i,k,l,m in zip(zlist,ngals[0,j,:],ngals[1,j,:],ngals[2,j,:],ngals[3,j,:],ngals[4,j,:],ngals[5,j,:],ngals[6,j,:],ngals[7,j,:],ngals[8,j,:],ngals[9,j,:],ngals[10,j,:]):
           print(a,b,c,d,e,f,g,h,i,k,l,m)
       ngals = nbrightgals
       j = 1
       print('# ')
       print('#redshift SMD1 SMD2 SMD3 SMD4 SMD5 SMD6 SMD7 SMD8 SMD9')
       for a,b,c,d,e,f,g,h,i,k in zip(zlist,ngals[0,j,:],ngals[1,j,:],ngals[2,j,:],ngals[3,j,:],ngals[4,j,:],ngals[5,j,:],ngals[6,j,:],ngals[7,j,:],ngals[8,j,:]):
           print(a,b,c,d,e,f,g,h,i,k)
       print('# ')
       print('#redshift SMD1(WE) SMD2(WE) SMD3(WE) SMD4(WE) SMD5(WE) SMD6(WE) SMD7(WE) SMD8(WE) SMD9(WE)')
       ngals = nbrightgals_err
       for a,b,c,d,e,f,g,h,i,k in zip(zlist,ngals[0,j,:],ngals[1,j,:],ngals[2,j,:],ngals[3,j,:],ngals[4,j,:],ngals[5,j,:],ngals[6,j,:],ngals[7,j,:],ngals[8,j,:]):
           print(a,b,c,d,e,f,g,h,i,k)
  
    plot_lir_lf(plt, outdir, obsdir, LFs_dust, file_name, csfr)
 
if __name__ == '__main__':
    main(*common.parse_args())
