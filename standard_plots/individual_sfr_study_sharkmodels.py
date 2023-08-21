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

import numpy as np
import utilities_statistics as us
import common

plt = common.load_matplotlib()

xtit="$\\rm LBT/Gyr$"
ytit="$\\rm log_{10}(SFR/M_{\odot} yr^{-1})$"

xmin, xmax, ymin, ymax = 0, 13.6, -1.5, 2.8
xleg = xmax + 0.025 * (xmax-xmin)
yleg = ymax - 0.07 * (ymax-ymin)

fig = plt.figure(figsize=(6.5,8))
mbins = ['9.125', '9.375', '9.625', '9.875', '10.125', '10.375', '10.625', '10.875', '11.125', '11.375', '11.625', '11.875']
colors = ('Navy','Blue','RoyalBlue','SkyBlue','Teal','DarkTurquoise','Aquamarine','Yellow', 'Gold',  'Orange','OrangeRed', 'LightSalmon', 'Crimson', 'Red', 'DarkRed')

ax = fig.add_subplot(111)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(2, 2, 0.2, 0.2))

for j in range(0,len(mbins)):
    print('../data/Models/SharkVariations/SFH/SFH_Lagos18_' + mbins[j] + '.dat')
    lbt_L18, sfr_L18 = np.loadtxt('../data/Models/SharkVariations/SFH/SFH_Lagos18_' + mbins[j] + '.dat', unpack = True, usecols = [0,1])
    lbt_L23, sfr_L23 = np.loadtxt('../data/Models/SharkVariations/SFH/SFH_Lagos23_' + mbins[j] + '.dat', unpack = True, usecols = [0,1])
    ax.plot(lbt_L23, sfr_L23, color=colors[j], linewidth=3, linestyle='solid', label='Shark v2.0'  if j == 0 else None) 
    ax.plot(lbt_L18, sfr_L18, color=colors[j], linewidth=3, linestyle='dotted', label = 'Shark v1.1 (L18)' if j == 0 else None)
    if(j < 3):
       print(j)
       ax.text(9.5,-1.1 + j * 0.2, mbins[j], fontsize=11, color=colors[j])
    if((j > 2) & (j<= 5)):
       print(j)
       ax.text(9.5,-0.4 + (j-3) * 0.3, mbins[j], fontsize=11, color=colors[j])
    if(j == 6):
       ax.text(9.5,0.55, mbins[j], fontsize=11, color=colors[j])
    if((j > 6) & (j <10)):
       print(j)
       ax.text(10.5,0.9 + (j-7) * 0.35, mbins[j], fontsize=11, color=colors[j])
    if(j==10):
       print(j)
       ax.text(10.5,2.15, mbins[j], fontsize=11, color=colors[j])


common.prepare_legend(ax, ['k','k'], loc=2)
plt.tight_layout()
common.savefig('/home/clagos/scm/git/shark/standard_plots', fig, "SFHs_BothSharkVersions.pdf")


### same plot but vs redshift:

xtit="$\\rm redshift$"
xmin, xmax, ymin, ymax = 0, 5, -1.5, 2.8
xleg = xmax + 0.025 * (xmax-xmin)
yleg = ymax - 0.07 * (ymax-ymin)

fig = plt.figure(figsize=(6.5,8))

ax = fig.add_subplot(111)
common.prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(2, 2, 0.2, 0.2))

for j in range(0,len(mbins)):
    print('../data/Models/SharkVariations/SFH/SFH_Lagos18_' + mbins[j] + '.dat')
    lbt_L18, sfr_L18 = np.loadtxt('../data/Models/SharkVariations/SFH/SFH_Lagos18_' + mbins[j] + '.dat', unpack = True, usecols = [0,1])
    lbt_L23, sfr_L23 = np.loadtxt('../data/Models/SharkVariations/SFH/SFH_Lagos23_' + mbins[j] + '.dat', unpack = True, usecols = [0,1])
    ax.plot(us.redshift(lbt_L23), sfr_L23, color=colors[j], linewidth=3, linestyle='solid', label='Shark v2.0'  if j == 0 else None) 
    ax.plot(us.redshift(lbt_L18), sfr_L18, color=colors[j], linewidth=3, linestyle='dotted', label = 'Shark v1.1 (L18)' if j == 0 else None)


common.prepare_legend(ax, ['k','k'], loc=2)
plt.tight_layout()
common.savefig('/home/clagos/scm/git/shark/standard_plots', fig, "SFHs_redshift_BothSharkVersions.pdf")

