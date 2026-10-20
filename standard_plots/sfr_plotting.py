
import math

import numpy as np

import common
import utilities_statistics as us

outdir = '/scratch/pawsey0119/clagos/SHARK_Out/Plots/medi-SURFS/Sharkv2-Lagos23/'

plt = common.load_matplotlib()

ztarget = 3.0
sfr_thresh = 10**(-10.75)
lbt_target = us.look_back_time(ztarget)


mgals, sfrgals, typeg = np.loadtxt('PassiveGalaxies_PMill.txt', unpack=True, usecols=[0,1,2])
ngals = len(sfrgals)

#read histories #########################
sfr_hist = np.loadtxt('SFHs_PassiveGalaxies_PMill.txt')
ms_hist = np.loadtxt('Masshistory_PassiveGalaxies_PMill.txt')

#select centrals only
cens = np.where((typeg == 0) & (sfrgals/mgals <= sfr_thresh))
ncens = len(typeg[cens])
mgals = mgals[cens]
sfrgals = sfrgals[cens]
typeg = typeg[cens]
ms_hist = ms_hist[cens,:]
sfr_hist = sfr_hist[cens,:]
ms_hist = ms_hist[0]
sfr_hist = sfr_hist[0]

np.savetxt('SFHs_PassiveCentrals_PMill.txt', sfr_hist)
np.savetxt('Masshistory_PassiveCentrals_PMill.txt', ms_hist)

props = np.zeros(shape = (ncens, 3))
props[:,0] = mgals
props[:,1] = sfrgals
props[:,2] = typeg
np.savetxt('PassiveCentrals_PMill.txt', props)

ids = np.argsort(mgals)
mgals = mgals[ids]
sfrgals = sfrgals[ids]
typeg = typeg[ids]
ms_hist = ms_hist[ids,:]
sfr_hist = sfr_hist[ids,:]

age_50 = np.zeros(shape = ncens)
age_80 = np.zeros(shape = ncens)
print("Number of galaxies:", ncens)

ind = np.where(sfr_hist == 0)
sfr_hist[ind] = 1e-4


#defineLBT bins ##############################################################################
#lbt_min = 11.7
#lbt_max = 13.7
#ntime = 50
#dtime = (lbt_max - lbt_min)/ntime
#tbins = np.arange(lbt_min, lbt_max, dtime)
#lbt_bins = tbins + dtime/2.0

lbt_bins = np.loadtxt('LBT_Lagos24.txt')
lbt_bins = lbt_bins - min(lbt_bins)
# measure stellar ages #######################################################################

def interpolate_ages(mtarget, m1, m2, lbt1, lbt2):
    m = (lbt2 - lbt1) / (m2-m1)
    y0 = lbt2 - m * m2
    return m * mtarget + y0


for g in range(0,ncens):
    for b in range(0,len(lbt_bins)):
        if((ms_hist[g,b] > 0.5 * mgals[g]) & (age_50[g] == 0)):
            age_50[g] = interpolate_ages ( 0.5 * mgals[g], ms_hist[g,b], ms_hist[g,b-1], lbt_bins[b], lbt_bins[b-1])
        if((ms_hist[g,b] > 0.8 * mgals[g]) & (age_80[g] == 0)):
            age_80[g] = interpolate_ages ( 0.8 * mgals[g], ms_hist[g,b], ms_hist[g,b-1], lbt_bins[b], lbt_bins[b-1])


print(age_50, age_80)
############## plot star formation histories ##################################################
xtit="$\\rm LBT/Gyr$"
ytit="$\\rm log_{10}(SFR/M_{\odot} yr^{-1})$"

xmin, xmax, ymin, ymax = 0, 2, -3.2, 3.5
xleg = xmax + 0.025 * (xmax-xmin)
yleg = ymax - 0.07 * (ymax-ymin)

fig = plt.figure(figsize=(6.5,5))
mbins = [10.0,10.1,10.2,10.3,10.4,10.5,10.75,12.5]
labels = ['10.05', '10.15', '10.25', '10.35', '10.45', '10.65', '11']
colors = ['Navy','SkyBlue','Aquamarine','Gold','Orange', 'LightSalmon','Red', 'DarkRed']

ax = fig.add_subplot(111)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.5, 0.5, 1, 1))

for j in range(0,len(mbins)-1):
    ax.text(0.05 + j*0.25, 3.6, labels[j], fontsize=11, color=colors[j])
    
for gal in range(0,ncens):
    select = np.where(mbins >= np.log10(mgals[gal]))
    indices = select[0]
    cols = colors[indices[0]-1]
    ax.plot(lbt_bins, np.log10(sfr_hist[gal,:]), color=cols, linewidth=1)

ax.plot([xmin,xmax], [np.log10(300.0),np.log10(300.0)], ls='solid',color='gray')
ax.text(xmin+0.1,np.log10(300.0)+0.1, "$\\rm 300\\, M_{\\odot}\\, yr^{-1}$", fontsize=11, color='gray')
ax.text(1.3,3, "SHARK", fontsize=13)
plt.tight_layout()
common.savefig(outdir, fig, "SFHs_massive_passive_galaxies_z3_centrals_Lagos24.pdf")


############### plot stellar mass growth histories #################################################################
ytit="$\\rm log_{10}(M_{\\rm stars}/M_{\odot})$"
ymin, ymax = 8,12
fig = plt.figure(figsize=(6.5,5))

ax = fig.add_subplot(111)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.5, 0.5, 1, 1))

for j in range(0,len(mbins)-1):
    ax.text(0.05 + j*0.25, 12.05, labels[j], fontsize=11, color=colors[j])

for gal in range(0,ncens):
    select = np.where(mbins >= np.log10(mgals[gal]))
    indices = select[0]
    cols = colors[indices[0]-1]
    ax.plot(lbt_bins, np.log10(ms_hist[gal,:]), color=cols, linewidth=1)
ax.text(1.3,11.7, "SHARK", fontsize=13)
ax.plot([xmin,xmax], [10,10], ls='solid',color='gray')

plt.tight_layout()
common.savefig(outdir, fig, "SMGrowth_massive_passive_galaxies_z3_centrals_Lagos24.pdf")



####################### plot age-mass relation ######################################################################
xtit="$\\rm log_{10}(M_{\\rm stars}/M_{\odot})$"
ytit="$\\rm age_{\\rm 50,90}/Gyr$"

xmin, xmax, ymin, ymax = 9.9, 12, 0, 2
xleg = xmax + 0.025 * (xmax-xmin)
yleg = ymax - 0.07 * (ymax-ymin)

fig = plt.figure(figsize=(6.5,5))
ax = fig.add_subplot(111)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.25, 0.25, 0.25, 0.25))

ax.plot(np.log10(mgals), age_50, ls=None, linewidth=0, marker='s', color='red', label = '$\\rm age_{50}$')
ax.plot(np.log10(mgals), age_80, ls=None, linewidth=0, marker='o', color='blue', label = '$\\rm age_{80}$')

m_in = np.log10(mgals)
for j in range(0,ncens):
    ax.plot([m_in[j],m_in[j]], [age_50[j],age_80[j]],ls='dotted', color='grey',linewidth=1)

ax.text(10,1.7, "SHARK", fontsize=13)
common.prepare_legend(ax, ['red','blue'], loc = 4)


plt.tight_layout()
common.savefig(outdir, fig, "age_mass_passive_galaxies_z3_centrals_Lagos24.pdf")

