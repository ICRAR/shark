import numpy as np
import scipy 
import utilities_statistics as us
import common

model = lagos13
redshift_power = 0.13405459588358698
v_sn = 120
beta_disk = 3.79746174188
eps_halo = 2.0
min_beta = 0.104050197191


vhot = v_sn * (age_univ)**redshift_power
const_sn =  (vhot/v)**power_index
age_univ =  us.look_back_time(redshifts)


zlow = 0
zupp = 20
dz = 0.125
zbins = np.arange(zlow,zupp,dz)
xzf = zbins + dz/2.0


plt = common.load_matplotlib()

fig = plt.figure(figsize=(5,4.5))
xtit = "$\\rm log_{10} (\\rm M_{\\star}/M_{\odot}) - 0.66 log_{10}(\\rm SFR/M_{\odot} yr^{-1})$"
ytit = "$\\rm log_{10}(\\rm Z_{\\rm gas}/Z_{\odot})$"
xmin, xmax, ymin, ymax = 7.1, 12, -3, 1
xleg = xmax - 0.2 * (xmax - xmin)
yleg = ymax - 0.1 * (ymax - ymin)

ax = fig.add_subplot(111)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(0.1, 1, 0.1, 1))
ax.text(xleg, yleg, 'z=0')

#Predicted relation
ind = np.where(fmzr[0,0,:] != 0)
yplot = (fmzr[0,0,ind])
errdn = (fmzr[0,1,ind])
errup = (fmzr[0,2,ind])

xplot = xmf[ind]
ax.errorbar(xplot,yplot[0],yerr=[errdn[0],errup[0]], ls='None', mfc='None', ecolor = 'b', mec='b',marker='o',label="all galaxies")

#observations approximate of Tremonti et al. (2004)
xG = [7.3,10.5]
yG = [-0.97,0.41]
yGu = [-0.77,0.61]
yGd = [-1.17,0.21]

ax.plot(xG,yG,'k',label="Andrews & Martini (2013)")
ax.plot(xG,yGu,'k',linestyle='dotted')
ax.plot(xG,yGd,'k',linestyle='dotted')

common.prepare_legend(ax, ['k','b'], loc=4)
common.savefig(outdir, fig, 'fmzr.pdf')

