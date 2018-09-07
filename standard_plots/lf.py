
#! /usr/bin/env python

# To use python modules in home directory
import h5py
import numpy as np
import os.path, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import scipy.optimize as so
import csv

import utilities_statistics as us 

plt.rcParams['legend.numpoints'] = 1

nom_lf_sdss_z0    = "SDSS_luminosity_functions.pdf"
nom_lf_ukidss_z0  = "UKIDSS_luminosity_functions.pdf"
nom_lf_IR_z0      = "IR_luminosity_functions.pdf"

# Initialize arguments
models= ["Shark-Lagos18"]
simu  = 'medi-SURFS'
imf   = 'cha'
dir = '/opt/claudia/SHArk_Out/' + simu + '/' 
dirphoto = '/Photometry/ProSpect/'
plotdir = '/opt/claudia/SHArk_Out/Plots/'
zlist = ["199"]

##################################
# Constants
GyrToYr = 1e9
Zsun = 0.0127
XH = 0.72
PI = 3.141592654
MpcToKpc = 1e3


##################################
# Mass function initialization
mlow = -30 + 5.0 * np.log10(0.677)
mupp = -10 + 5.0 * np.log10(0.677)
dm = 0.5
mbins = np.arange(mlow,mupp,dm)
xlf   = mbins + dm/2.0 

# colour distribution initialization
clow  = -0.1
cupp  = 1.5
dc    = 0.5
cbins = np.arange(clow,cupp,dc)
xc    = cbins + dc/2.0 

magbins = [-17.13,-17.88,-18.63,-19.38,-20.13,-20.88,-21.63]

hist_ur       = np.zeros(shape = (len(magbins)-1, len(cbins)))
hist_gr       = np.zeros(shape = (len(magbins)-1, len(cbins)))


#Flux corrections for dust from Hill+2010 for bands ugrizYJHK
flux_corrections = [2.27, 1.69, 1.56, 1.47, 1.41, 1.35, 1.32, 1.22, 1.15]
mag_corrections  = -2.5*np.log10(flux_corrections)
#structure of CSV:
#bands = [FUV, NUV, uS, gS, rS, iS, ZV, YV, JV, HV, KV, W1, W2, W3, W4, P100, P160, S250, S350, S500] for bulges built via disk ins, mergers, total bulges, disks and total galaxy. (20 bands) 
#ID,
#bulge formed via disk instabilities no dust absolute magnitude
#ab_mag_nodust
#ap_mag_nodust
#ab_mag_dust
#ap_mag_dust
# Loop over redshift and subvolumes
for model in range(0,len(models)):
	for index in range(0,len(zlist)):
		indir = dir + models[model] + '/' + dirphoto + zlist[index] + '/' + '0/'
		file = indir+'SharkSED.csv'
		with open(file) as csvfile:
			reader = csv.reader(csvfile, delimiter=',')
			ngal = 0L
			nbands = 0
			for row in reader:
				#define number of bands only once (based on a number of entries of 3 galaxy components x 4 mag/app/dust/no-dust, which gives 12).
				if(ngal == 0):
					nbands = (len(row)-1) / (5*2*2)
				#count number of galaxies
				ngal = ngal + 1
			ngal = ngal -1 #ignore first row
			#allocate SED arrays.
			print 'number of galaxies', ngal,' number of bands',nbands
			IDs       = np.zeros(shape = (ngal))
			ID_Shark  = np.zeros(shape = (ngal))
			SED 	  = np.zeros(shape = (ngal, nbands, 5, 2, 2)) #5: bulge disk-instabilities, bulge mergers, bulge, disk and total; 2:absolute and apparent magnitude; 2: no dust and dust.

                with open(file) as csvfile:
			j=0L
			lines = 0
			reader = csv.reader(csvfile, delimiter=',')
		
			for row in reader:
				if(lines > 0):
					c	 =  0L
					for d in range(0,2):
						for mag in range (0,2):
							for comp in range(0,5):
								for band in range (0, nbands):
									IDs[j]      = row[0]
                                                                        ID_Shark[j] = row[1]
									magin       = row[c+2]
									if(magin):
										SED[j,band,comp,mag,d] = float(magin)
									c = c+1
					j = j+1		
				lines = lines +1

		if(model == 0):
			hist_lf            = np.zeros(shape = (len(models), len(zlist), len(mbins), nbands))
			hist_lf_disk       = np.zeros(shape = (len(models), len(zlist), len(mbins), nbands))
			hist_lf_bulge      = np.zeros(shape = (len(models), len(zlist), len(mbins), nbands))
			hist_lf_pseudobulge= np.zeros(shape = (len(models), len(zlist), len(mbins), nbands))

			hist_lf_nodust            = np.zeros(shape = (len(models), len(zlist), len(mbins), nbands))
			hist_lf_disk_nodust       = np.zeros(shape = (len(models), len(zlist), len(mbins), nbands))
			hist_lf_bulge_nodust      = np.zeros(shape = (len(models), len(zlist), len(mbins), nbands))
			hist_lf_pseudobulge_nodust= np.zeros(shape = (len(models), len(zlist), len(mbins), nbands))

			#pow(10.0,mag/(-2.5))
		for i in range(0,nbands):
			#dust magnitude total
			SED[:,i,4,0,1] = np.log10(pow(10.0,SED[:,i,2,0,1]/(-2.5)) + pow(10.0,SED[:,i,3,0,1]/(-2.5)) ) * (-2.5)
			SED[:,i,4,0,0] = np.log10(pow(10.0,SED[:,i,2,0,0]/(-2.5)) + pow(10.0,SED[:,i,3,0,0]/(-2.5)) ) * (-2.5)

	                ind = np.where(SED[:,i,4,0,1] < -1)
	                H, bins_edges = np.histogram(SED[ind,i,4,0,1],bins=np.append(mbins,mupp))
	                hist_lf[model,index,:,i] = hist_lf[model,index,:,i] + H
	                ind = np.where(SED[:,i,0,0,1] < -1)
	                H, bins_edges = np.histogram(SED[ind,i,0,0,1],bins=np.append(mbins,mupp))
	                hist_lf_pseudobulge[model,index,:,i] = hist_lf_pseudobulge[model,index,:,i] + H
	                ind = np.where(SED[:,i,1,0,1] < -1)
	                H, bins_edges = np.histogram(SED[ind,i,1,0,1],bins=np.append(mbins,mupp))
	                hist_lf_bulge[model,index,:,i] = hist_lf_bulge[model,index,:,i] + H
	                ind = np.where(SED[:,i,3,0,1] < -1)
	                H, bins_edges = np.histogram(SED[ind,i,3,0,1],bins=np.append(mbins,mupp))
	                hist_lf_disk[model,index,:,i] = hist_lf_disk[model,index,:,i] + H


			#no dust magnitudes
	                ind = np.where(SED[:,i,4,0,0] < -1)
	                H, bins_edges = np.histogram(SED[ind,i,4,0,0],bins=np.append(mbins,mupp))
	                hist_lf_nodust[model,index,:,i] = hist_lf_nodust[model,index,:,i] + H
	                ind = np.where(SED[:,i,1,0,0] < -1)
	                H, bins_edges = np.histogram(SED[ind,i,1,0,0],bins=np.append(mbins,mupp))
	                hist_lf_bulge_nodust[model,index,:,i] = hist_lf_bulge_nodust[model,index,:,i] + H
	                ind = np.where(SED[:,i,0,0,0] < -1)
	                H, bins_edges = np.histogram(SED[ind,i,0,0,0],bins=np.append(mbins,mupp))
	                hist_lf_pseudobulge_nodust[model,index,:,i] = hist_lf_pseudobulge_nodust[model,index,:,i] + H
	                ind = np.where(SED[:,i,3,0,0] < -1)
	                H, bins_edges = np.histogram(SED[ind,i,3,0,0],bins=np.append(mbins,mupp))
	                hist_lf_disk_nodust[model,index,:,i] = hist_lf_disk_nodust[model,index,:,i] + H

		#read volume and h factor
                indir = dir + models[model] + '/' + zlist[index] + '/' + '0/'
                file = indir+'galaxies.hdf5'

                if (os.path.isfile(file)):
                        f = h5py.File(file, 'r')
                        group = f['Galaxies']
                        h0 = f['Cosmology/h'].value
                        volh = f['runInfo/EffectiveVolume'].value
                print volh 

for model in range(0,len(models)):
        for index in range(0,len(zlist)):
		for band in range(0,nbands):
	                if(volh > 0.):
        	                vol = volh  # In Mpc^3
                	        hist_lf[model,index,:,band]        = hist_lf[model,index,:,band]/vol
                        	hist_lf_disk[model,index,:,band]   = hist_lf_disk[model,index,:,band]/vol
	                        hist_lf_bulge[model,index,:,band]  = hist_lf_bulge[model,index,:,band]/vol
	                        hist_lf_pseudobulge[model,index,:,band]  = hist_lf_pseudobulge[model,index,:,band]/vol

                	        hist_lf_nodust[model,index,:,band]        = hist_lf_nodust[model,index,:,band]/vol
                        	hist_lf_disk_nodust[model,index,:,band]   = hist_lf_disk_nodust[model,index,:,band]/vol
	                        hist_lf_bulge_nodust[model,index,:,band]  = hist_lf_bulge_nodust[model,index,:,band]/vol
	                        hist_lf_pseudobulge_nodust[model,index,:,band]  = hist_lf_pseudobulge_nodust[model,index,:,band]/vol

                        	# Take logs
	                        ind = np.where(hist_lf[model,index,:,band] > 0.)
	                        hist_lf[model,index,ind,band]  =  np.log10(hist_lf[model,index,ind,band])
        	        	ind = np.where(hist_lf_bulge[model,index,:,band] > 0.)
               	       		hist_lf_bulge[model,index,ind,band]  =  np.log10(hist_lf_bulge[model,index,ind,band])
                        	ind = np.where(hist_lf_disk[model,index,:,band] > 0.)
                        	hist_lf_disk[model,index,ind,band] = np.log10(hist_lf_disk[model,index,ind,band])
         	        	ind = np.where(hist_lf_pseudobulge[model,index,:,band] > 0.)
               	       		hist_lf_pseudobulge[model,index,ind,band]  =  np.log10(hist_lf_pseudobulge[model,index,ind,band])

	                        ind = np.where(hist_lf_nodust[model,index,:,band] > 0.)
	                        hist_lf_nodust[model,index,ind,band]  =  np.log10(hist_lf_nodust[model,index,ind,band])
        	        	ind = np.where(hist_lf_bulge_nodust[model,index,:,band] > 0.)
               	       		hist_lf_bulge_nodust[model,index,ind,band]  =  np.log10(hist_lf_bulge_nodust[model,index,ind,band])
                        	ind = np.where(hist_lf_disk_nodust[model,index,:,band] > 0.)
                        	hist_lf_disk_nodust[model,index,ind,band] = np.log10(hist_lf_disk_nodust[model,index,ind,band])
        	        	ind = np.where(hist_lf_pseudobulge_nodust[model,index,:,band] > 0.)
               	       		hist_lf_pseudobulge_nodust[model,index,ind,band]  =  np.log10(hist_lf_pseudobulge_nodust[model,index,ind,band])


#compute colour distributions
uband = 2
gband = 3
rband = 4
ubandl = SED[ind,uband,4,0,1]
gbandl = SED[ind,gband,4,0,1]
rbandl = SED[ind,rband,4,0,1]

#hist_ur       = np.zeros(shape = (len(magbins)-1, len(cbins)))
#hist_gr       = np.zeros(shape = (len(magbins)-1, len(cbins)))

for mag in range(0,len(magbins)-1):
    ind = np.where((ubandl < -1) && (gbandl < -1) && (rbandl < magbins[mag]) && (bandl >= magbins[mag+1]))
    H, bins_edges  = np.histogram(ubandl[ind] - rbandl[ind],bins=np.append(cbins,cupp))
    hist_ur[mag,:] = hist_ur[mag,:] + H
    hist_ur[mag,:] = hist_ur[mag,:] / (len(ubandl[ind]) * dc)
    H, bins_edges  = np.histogram(gbandl[ind] - rbandl[ind],bins=np.append(cbins,cupp))
    hist_gr[mag,:] = hist_gr[mag,:] + H
    hist_gr[mag,:] = hist_gr[mag,:] / (len(gbandl[ind]) * dc)


###################################
#Observations directory
obsdir = '/opt/claudia/ObsData/lf/'


#plot colour distributions
fig = plt.figure(figsize=(9.7,11.7))
xtit = "$(u-r)$"
ytit = "$\\rm dp/d(u-r)$"
xmin, xmax, ymin, ymax = 0, 3.5, 0, 2.5
xleg = xmax - 0.2 * (xmax - xmin)
yleg = ymax - 0.1 * (ymax - ymin)

subplots = (321, 322, 323, 324, 325, 326)
indeces = (0, 1, 2, 3, 4, 5)
for subplot, idx in zip(subplots, indeces):

    ax = fig.add_subplot(subplot)
    common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1))

    # Observations
    for obs, marker in obs_and_markers:
        common.errorbars(ax, obs.x, obs.y, obs.yerrdn, obs.yerrup, 'grey',
                         marker, err_absolute=obs.err_absolute, label=obs.label)

    # Predicted SMF
    if plotz[idx]:
        y = hist_smf[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'r', label='all galaxies' if idx == 0 else None)
        y = hist_smf_cen[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'b', linestyle='dotted', label ='centrals' if idx == 0 else None)
        y = hist_smf_sat[idx,:]
        ind = np.where(y < 0.)
        ax.plot(xmf[ind],y[ind],'g', linestyle='dashed', label ='satellites' if idx == 0 else None)

        if z < 1:
            y = hist_smf_30kpc[idx,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'k', linestyle='dotted', linewidth=1, label ='30kpc')
        if z >= 1:
            y = hist_smf_err[idx,:]
            ind = np.where(y < 0.)
            ax.plot(xmf[ind],y[ind],'r', linestyle='dashdot', linewidth=2, label ='0.25 dex error')

    colors = []
    if idx == 0:
        colors = ['r','b','g']
    if z < 1:
        colors += ['k']
    if z >= 1:
        colors = ['r']
    colors += ['grey', 'grey','grey']

    common.prepare_legend(ax, colors)

common.savefig(outdir, fig, 'stellarmf_z.pdf')


###################################
#   Plots stellar mass frunction
fig = plt.figure(figsize=(12,12))

#fig.suptitle(model+' IMF='+imf,fontsize = 9)

xtit="$\\rm mag-5log(h) (AB)$"
ytit="$\\rm log_{10}(\Phi/(0.5\\, {\\rm mag})/h^3 {\\rm Mpc}^{-3})$"

xmin = -25
xmax = -13
ymin = -5.1
ymax = -1.

# Legend
xleg= xmin + 0.2*(xmax-xmin)
yleg= ymax - 0.1*(ymax-ymin)

fs = 13
hcorr = 5.0*np.log10(0.677)
xlf_obs  = xlf - hcorr 

# z=0 ##################################
iz=0
ax = fig.add_subplot(331)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.set_ylabel(ytit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 GALEX FUV-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
file = obsdir+'lf1500_z0_driver12.data'
lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
indx = np.where(p > 0)
yobs = np.log10(p[indx])
ydn  = np.log10(p[indx]-dp[indx])
yup  = np.log10(p[indx]+dp[indx])

ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Driver+2012")

###################################
###################################
###################################
###################################
#Predicted LF

band = 0
ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

###################################
ax = fig.add_subplot(332)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 GALEX NUV-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
file = obsdir+'lf2300_z0_driver12.data'
lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
indx = np.where(p > 0)
yobs = np.log10(p[indx])
ydn  = np.log10(p[indx]-dp[indx])
yup  = np.log10(p[indx]+dp[indx])

ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

###################################
###################################
###################################
###################################
#Predicted LF
band = 1

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')


# z=0 ##################################
iz=0
ax = fig.add_subplot(334)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.set_ylabel(ytit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 SDSS u-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
file = obsdir+'lfu_z0_driver12.data'
lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
indx = np.where(p > 0)
yobs = np.log10(p[indx])
ydn  = np.log10(p[indx]-dp[indx])
yup  = np.log10(p[indx]+dp[indx])

ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Driver+2012")

###################################
###################################
###################################
###################################
#Predicted LF

band = 2
ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3, label ='Shark')
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1, label ='Shark intrinsic')

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted', label ='disks')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed', label ='bulges')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot',label='disk-ins-driven bulges')

# Legend
leg = ax.legend(bbox_to_anchor=[3.2,-0.3],prop={'size':12})
colors = ['k','k','b','r','LightSalmon','grey','grey']
for color,text in zip(colors,leg.get_texts()):
    text.set_color(color)
leg.draw_frame(False)

###################################
ax = fig.add_subplot(335)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 SDSS g-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
file = obsdir+'lfg_z0_driver12.data'
lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
indx = np.where(p > 0)
yobs = np.log10(p[indx])
ydn  = np.log10(p[indx]-dp[indx])
yup  = np.log10(p[indx]+dp[indx])

ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

###################################
###################################
###################################
###################################
#Predicted LF
band = 3 

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

###################################
ax = fig.add_subplot(336)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.set_xlabel(xtit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 SDSS r-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
file = obsdir+'lfr_z0_driver12.data'
lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
indx = np.where(p > 0)
yobs = np.log10(p[indx])
ydn  = np.log10(p[indx]-dp[indx])
yup  = np.log10(p[indx]+dp[indx])

ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

###################################
###################################
###################################
###################################
#Predicted LF
band = 4 

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

###################################
ax = fig.add_subplot(337)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.set_xlabel(xtit,fontsize = fs)
ax.set_ylabel(ytit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 SDSS i-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
file = obsdir+'lfi_z0_driver12.data'
lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
indx = np.where(p > 0)
yobs = np.log10(p[indx])
ydn  = np.log10(p[indx]-dp[indx])
yup  = np.log10(p[indx]+dp[indx])

ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

###################################
###################################
###################################
###################################
#Predicted LF
band = 5

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

###################################
ax = fig.add_subplot(338)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.set_xlabel(xtit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 SDSS z-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
file = obsdir+'lfz_z0_driver12.data'
lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
indx = np.where(p > 0)
yobs = np.log10(p[indx])
ydn  = np.log10(p[indx]-dp[indx])
yup  = np.log10(p[indx]+dp[indx])

ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

###################################
###################################
###################################
###################################
#Predicted LF
band = 6

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')


###################################
#Save figure
print 'making plot '+ nom_lf_sdss_z0
plotfile = plotdir + simu + '/' + models[0] + '/' + nom_lf_sdss_z0 
fig.savefig(plotfile, dvi=300, pad_inches=0)

##################################
##################################
#UKIDSS
###################################
#   Plots stellar mass frunction
fig = plt.figure(figsize=(12,8.5))

#fig.suptitle(model+' IMF='+imf,fontsize = 9)

xtit="$\\rm mag-5log(h) (AB)$"
ytit="$\\rm log_{10}(\Phi/(0.5\\, {\\rm mag})/h^3 {\\rm Mpc}^{-3})$"

xmin = -25
xmax = -14
ymin = -5.1
ymax = -1.

# Legend
xleg= xmin + 0.2*(xmax-xmin)
yleg= ymax - 0.1*(ymax-ymin)

fs = 13

# z=0 ##################################
iz=0
ax = fig.add_subplot(231)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.set_ylabel(ytit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 UKIDSS Y-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
file = obsdir+'lfy_z0_driver12.data'
lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
indx = np.where(p > 0)
yobs = np.log10(p[indx])
ydn  = np.log10(p[indx]-dp[indx])
yup  = np.log10(p[indx]+dp[indx])

ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Driver+2012")

###################################
###################################
###################################
###################################
#Predicted LF

band = 7
ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3, label ='Shark')
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1,label ='Shark intrinsic')

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted', label ='disks')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed', label ='merger-driven bulges')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot', label='disk-ins-driven bulges')

# Legend
leg = ax.legend(loc=4,prop={'size':12})
leg = ax.legend(bbox_to_anchor=[3.2,-0.3],prop={'size':12})
colors = ['k','k','b','r','LightSalmon','grey','grey']
for color,text in zip(colors,leg.get_texts()):
    text.set_color(color)
leg.draw_frame(False)

###################################
ax = fig.add_subplot(232)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 UKIDSS J-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
file = obsdir+'lfj_z0_driver12.data'
lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
indx = np.where(p > 0)
yobs = np.log10(p[indx])
ydn  = np.log10(p[indx]-dp[indx])
yup  = np.log10(p[indx]+dp[indx])

ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

###################################
###################################
###################################
###################################
#Predicted LF
band = 8 

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

###################################
ax = fig.add_subplot(234)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.set_xlabel(xtit,fontsize = fs)
ax.set_ylabel(ytit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 UKIDSS H-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
file = obsdir+'lfh_z0_driver12.data'
lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
indx = np.where(p > 0)
yobs = np.log10(p[indx])
ydn  = np.log10(p[indx]-dp[indx])
yup  = np.log10(p[indx]+dp[indx])

ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

###################################
###################################
###################################
###################################
#Predicted LF
band = 9

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

###################################
ax = fig.add_subplot(235)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(xtit,fontsize = fs)
#ax.set_ylabel(ytit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 UKIDSS K-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
file = obsdir+'lfk_z0_driver12.data'
lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
indx = np.where(p > 0)
yobs = np.log10(p[indx])
ydn  = np.log10(p[indx]-dp[indx])
yup  = np.log10(p[indx]+dp[indx])

ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

###################################
###################################
###################################
###################################
#Predicted LF
band = 10

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')


###################################
#Save figure
print 'making plot '+ nom_lf_ukidss_z0
plotfile = plotdir + simu + '/' + models[0] + '/' + nom_lf_ukidss_z0
fig.savefig(plotfile, dvi=300, pad_inches=0)

##################################
##################################
#IR
###################################
#   Plots stellar mass frunction
fig = plt.figure(figsize=(12,12))

#fig.suptitle(model+' IMF='+imf,fontsize = 9)
#W2, W3, W4, P100, P160, S250, S350, S500

xtit="$\\rm mag-5log(h) (AB)$"
ytit="$\\rm log_{10}(\Phi/dex^{-1} h^3 {\\rm Mpc}^{-3})$"

xmin = -30
xmax = -16
ymin = -5.1
ymax = -1.

# Legend
xleg= xmin + 0.2*(xmax-xmin)
yleg= ymax - 0.1*(ymax-ymin)

fs = 13

# z=0 ##################################
iz=0
ax = fig.add_subplot(331)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.set_ylabel(ytit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 WISE 2-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
#file = obsdir+'lfy_z0_driver12.data'
#lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
#indx = np.where(p > 0)
#yobs = np.log10(p[indx])
#ydn  = np.log10(p[indx]-dp[indx])
#yup  = np.log10(p[indx]+dp[indx])

#ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Driver+2012")

###################################
###################################
###################################
###################################
#Predicted LF

band = 12
ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band] - np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band] - np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

###################################
ax = fig.add_subplot(332)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 WISE 3-band')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
#file = obsdir+'lfj_z0_driver12.data'
#lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
#indx = np.where(p > 0)
#yobs = np.log10(p[indx])
#ydn  = np.log10(p[indx]-dp[indx])
#yup  = np.log10(p[indx]+dp[indx])
#
#ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

###################################
###################################
###################################
###################################
#Predicted LF
band = 13

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band] - np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

###################################
ax = fig.add_subplot(333)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
#ax.set_xlabel(xtit,fontsize = fs)
#ax.set_ylabel(ytit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 WISE 4-band')

###################################
###################################
###################################
###################################
#Predicted LF
band = 14

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band] - np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')


file = obsdir+'L24microns.dat'
lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
dml       = 0.15
xm        = [25.0,25.0]
xm[1]     = xm[1] + dml
tenpctocm = 3.086e19
corrpc    = np.log10(tenpctocm)*2.0
xobs      = -2.5*(lm[0:13]-corrpc+7.0-np.log10(4.0*PI)+ np.log10(24.0/22.0) ) - 48.6 - hcorr 
xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)+ np.log10(24.0/22.0) ) - 48.6 - hcorr
dmm       = np.abs(xm[0] - xm[1])
yobs      = np.log10(pow(10.0,p[0:13]) * dml/dmm)-3.0*np.log10(0.677)
ydn       = np.log10(pow(10.0,dpn[0:13]) * dml/dmm)-3.0*np.log10(0.677)
yup       = np.log10(pow(10.0,dpu[0:13]) * dml/dmm)-3.0*np.log10(0.677)

#ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Marchetti+2016")

dml       = 0.4
xm        = [25.0,25.0]
xm[1]     = xm[1] + dml
tenpctocm = 3.086e19
corrpc    = np.log10(tenpctocm)*2.0
xobs      = -2.5*(lm[14:22]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr 
xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
dmm       = np.abs(xm[0] - xm[1])
yobs      = np.log10(pow(10.0,p[14:22]) * dml/dmm)-3.0*np.log10(0.677)
ydn       = np.log10(pow(10.0,dpn[14:22]) * dml/dmm)-3.0*np.log10(0.677)
yup       = np.log10(pow(10.0,dpu[14:22]) * dml/dmm)-3.0*np.log10(0.677)

#ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Marleau+2007")

dml       = 0.4
xm        = [25.0,25.0]
xm[1]     = xm[1] + dml
tenpctocm = 3.086e19
corrpc    = np.log10(tenpctocm)*2.0
xobs      = -2.5*(lm[23:28]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr 
xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
dmm       = np.abs(xm[0] - xm[1])
yobs      = np.log10(pow(10.0,p[23:28]) * dml/dmm)-3.0*np.log10(0.677)
ydn       = np.log10(pow(10.0,dpn[23:28]) * dml/dmm)-3.0*np.log10(0.677)
yup       = np.log10(pow(10.0,dpu[23:28]) * dml/dmm)-3.0*np.log10(0.677)

#ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Rodighiero+2010")

#leg = ax.legend(loc=4,prop={'size':12})
#colors = ['grey','grey','grey']
#for color,text in zip(colors,leg.get_texts()):
#    text.set_color(color)
#leg.draw_frame(False)


###################################
ax = fig.add_subplot(334)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.set_ylabel(ytit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 P100')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
#file = obsdir+'lfk_z0_driver12.data'
#lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
#indx = np.where(p > 0)
#yobs = np.log10(p[indx])
#ydn  = np.log10(p[indx]-dp[indx])
#yup  = np.log10(p[indx]+dp[indx])
#
#ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

###################################
###################################
###################################
###################################
#Predicted LF
band = 15

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band] - np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
y[:] = 0.0
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')



###################################
ax = fig.add_subplot(335)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 P160')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
#file = obsdir+'lfk_z0_driver12.data'
#lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
#indx = np.where(p > 0)
#yobs = np.log10(p[indx])
#ydn  = np.log10(p[indx]-dp[indx])
#yup  = np.log10(p[indx]+dp[indx])
#
#ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

###################################
###################################
###################################
###################################
#Predicted LF
band = 16

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')


file = obsdir+'L160microns.dat'
lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
dml       = 0.15
xm        = [25.0,25.0]
xm[1]     = xm[1] + dml
tenpctocm = 3.086e19
corrpc    = np.log10(tenpctocm)*2.0
xobs      = -2.5*(lm[0:14]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr 
xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
dmm       = np.abs(xm[0] - xm[1])
yobs      = np.log10(pow(10.0,p[0:14]) * dml/dmm)-3.0*np.log10(0.677)
ydn       = np.log10(pow(10.0,dpn[0:14]) * dml/dmm)-3.0*np.log10(0.677)
yup       = np.log10(pow(10.0,dpu[0:14]) * dml/dmm)-3.0*np.log10(0.677)

ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*')

dml       = 0.4
xm        = [25.0,25.0]
xm[1]     = xm[1] + dml
tenpctocm = 3.086e19
corrpc    = np.log10(tenpctocm)*2.0
xobs      = -2.5*(lm[15:21]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr 
xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
dmm       = np.abs(xm[0] - xm[1])
yobs      = np.log10(pow(10.0,p[15:21]) * dml/dmm)-3.0*np.log10(0.677)
ydn       = np.log10(pow(10.0,dpn[15:21]) * dml/dmm)-3.0*np.log10(0.677)
yup       = np.log10(pow(10.0,dpu[15:21]) * dml/dmm)-3.0*np.log10(0.677)

ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='p')


###################################
ax = fig.add_subplot(336)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
#ax.set_ylabel(ytit,fontsize = fs)
ax.set_xlabel(xtit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 S250')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
file = obsdir+'lf250_dye10.data'
lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
lx        = np.log10(lm[0:7])
dml       = np.abs(lx[1] - lx[0])
tenpctocm = 3.086e19
corrpc    = np.log10(tenpctocm)*2.0
xobs      = -2.5*(np.log10(lm[0:7])-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr 
dmm       = np.abs(xobs[0] - xobs[1])
yobs      = np.log10(p[0:7] * dml/dmm)-3.0*np.log10(0.677)
ydn       = np.log10((p[0:7]- dp[0:7]) * dml/dmm)-3.0*np.log10(0.677)
yup       = np.log10((p[0:7]+ dp[0:7]) * dml/dmm)-3.0*np.log10(0.677)

ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

file = obsdir+'L250microns.dat'
lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
dml       = 0.15
xm        = [25.0,25.0]
xm[1]     = xm[1] + dml
tenpctocm = 3.086e19
corrpc    = np.log10(tenpctocm)*2.0
xobs      = -2.5*(lm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr 
pxm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
dmm       = np.abs(xm[0] - xm[1])
yobs      = np.log10(pow(10.0,p) * dml/dmm)-3.0*np.log10(0.677)
ydn       = np.log10(pow(10.0,dpn) * dml/dmm)-3.0*np.log10(0.677)
yup       = np.log10(pow(10.0,dpu) * dml/dmm)-3.0*np.log10(0.677)

ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*')


###################################
###################################
###################################
###################################
#Predicted LF
band = 17

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

leg = ax.legend(loc=4,prop={'size':12})
colors = ['grey','grey']
for color,text in zip(colors,leg.get_texts()):
    text.set_color(color)
leg.draw_frame(False)


###################################
ax = fig.add_subplot(337)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.set_xlabel(xtit,fontsize = fs)
ax.set_ylabel(ytit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 S350')

###################################
###################################
###################################
###################################
#Predicted LF
band = 18

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3)

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot')

file = obsdir+'L350microns.dat'
lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
dml       = 0.2
xm        = [25.0,25.0]
xm[1]     = xm[1] + dml
tenpctocm = 3.086e19
corrpc    = np.log10(tenpctocm)*2.0
xobs      = -2.5*(lm[0:9]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr 
xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
dmm       = np.abs(xm[0] - xm[1])
yobs      = np.log10(pow(10.0,p[0:9]) * dml/dmm)-3.0*np.log10(0.677)
ydn       = np.log10(pow(10.0,dpn[0:9]) * dml/dmm)-3.0*np.log10(0.677)
yup       = np.log10(pow(10.0,dpu[0:9]) * dml/dmm)-3.0*np.log10(0.677)

ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*')

dml       = 0.3
xm        = [25.0,25.0]
xm[1]     = xm[1] + dml
tenpctocm = 3.086e19
corrpc    = np.log10(tenpctocm)*2.0
xobs      = -2.5*(lm[10:19]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr 
xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
dmm       = np.abs(xm[0] - xm[1])
yobs      = np.log10(pow(10.0,p[10:19]) * dml/dmm)-3.0*np.log10(0.677)
ydn       = np.log10(pow(10.0,dpn[10:19]) * dml/dmm)-3.0*np.log10(0.677)
yup       = np.log10(pow(10.0,dpu[10:19]) * dml/dmm)-3.0*np.log10(0.677)

ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s')


###################################
ax = fig.add_subplot(338)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax) 
ax.set_xlabel(xtit,fontsize = fs)
ax.xaxis.set_major_locator(MultipleLocator(2.)) 
ax.xaxis.set_minor_locator(MultipleLocator(2)) 
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1.)) 
ax.tick_params(labelsize=13)
ax.text(xleg,yleg, 'z=0 S500')

#Baldry (Chabrier IMF), ['Baldry+2012, z<0.06']
#file = obsdir+'lfk_z0_driver12.data'
#lm,p,dp = np.loadtxt(file,usecols=[0,1,2],unpack=True)
#indx = np.where(p > 0)
#yobs = np.log10(p[indx])
#ydn  = np.log10(p[indx]-dp[indx])
#yup  = np.log10(p[indx]+dp[indx])
#
#ax.errorbar(lm[indx], yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o')

###################################
###################################
###################################
###################################
#Predicted LF
band = 19

ind = np.where(hist_lf[0,0,:,band] < 0.)
y = hist_lf[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'k', linewidth=3, label ='Shark')
ind = np.where(hist_lf_nodust[0,0,:,band] < 0.)
y = hist_lf_nodust[0,0,ind,band]
y[:] = 0.0
ax.plot(xlf_obs[ind],y[0],'k', linewidth=1, label='Shark intrinsic')

ind = np.where(hist_lf_disk[0,0,:,band] < 0.)
y = hist_lf_disk[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'b', linewidth=2, linestyle='dotted', label ='disks')
ind = np.where(hist_lf_bulge[0,0,:,band] < 0.)
y = hist_lf_bulge[0,0,ind,band]- np.log10(dm)
ax.plot(xlf_obs[ind],y[0],'r', linewidth=2, linestyle='dashed',  label ='merger-driven bulges')
ind = np.where(hist_lf_pseudobulge[0,0,:,band] < 0.)
y = hist_lf_pseudobulge[0,0,ind,band]
ax.plot(xlf_obs[ind],y[0],'LightSalmon', linewidth=2, linestyle='dashdot', label='disk-ins-driven bulges')

file = obsdir+'L500microns.dat'
lm,p,dpn,dpu = np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
dml       = 0.2
xm        = [25.0,25.0]
xm[1]     = xm[1] + dml
tenpctocm = 3.086e19
corrpc    = np.log10(tenpctocm)*2.0
xobs      = -2.5*(lm[0:11]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr 
xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
dmm       = np.abs(xm[0] - xm[1])
yobs      = np.log10(pow(10.0,p[0:11]) * dml/dmm)-3.0*np.log10(0.677)
ydn       = np.log10(pow(10.0,dpn[0:11]) * dml/dmm)-3.0*np.log10(0.677)
yup       = np.log10(pow(10.0,dpu[0:11]) * dml/dmm)-3.0*np.log10(0.677)

ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='*',label="Marchetti+2016")

dml       = 0.3
xm        = [25.0,25.0]
xm[1]     = xm[1] + dml
tenpctocm = 3.086e19
corrpc    = np.log10(tenpctocm)*2.0
xobs      = -2.5*(lm[11:19]-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr 
xm        = -2.5*(xm-corrpc+7.0-np.log10(4.0*PI)) - 48.6 - hcorr
dmm       = np.abs(xm[0] - xm[1])
yobs      = np.log10(pow(10.0,p[11:19]) * dml/dmm)-3.0*np.log10(0.677)
ydn       = np.log10(pow(10.0,dpn[11:19]) * dml/dmm)-3.0*np.log10(0.677)
yup       = np.log10(pow(10.0,dpu[11:19]) * dml/dmm)-3.0*np.log10(0.677)

ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='s',label="Negrello+2013")

xobs = np.zeros(shape = 2)
xobs[0] = 1.0
xobs[1] = 1.0
yobs=xobs
ydn=xobs
yup=xobs
ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='p',label="Patel+2013")
ax.errorbar(xobs, yobs, yerr=[yobs-ydn,yup-yobs], ls='None', mfc='None', ecolor = 'grey', mec='grey',marker='o',label="Dye+2010")

# Legend
leg = ax.legend(bbox_to_anchor=[1.1,0.8],prop={'size':12})
colors = ['k','k','b','r','LightSalmon','grey','grey','grey','grey']
for color,text in zip(colors,leg.get_texts()):
    text.set_color(color)
leg.draw_frame(False)


###################################
#Save figure
print 'making plot '+ nom_lf_IR_z0
plotfile = plotdir + simu + '/' + models[0] + '/' + nom_lf_IR_z0
fig.savefig(plotfile, dvi=300, pad_inches=0)


