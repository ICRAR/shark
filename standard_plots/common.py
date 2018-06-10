#
#    ICRAR - International Centre for Radio Astronomy Research
#    (c) UWA - The University of Western Australia, 2018
#    Copyright by UWA (in the framework of the ICRAR)
#    All rights reserved
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#    MA 02111-1307  USA
#
"""Common routines for shark plots"""

import argparse
import collections
import os
import sys

import h5py
import numpy as np

PY2 = sys.version_info[0] == 2
if PY2:
    import ConfigParser as configparser
else:
    import configparser


def load_matplotlib():
    import matplotlib.pyplot as plt
    plt.rcParams['legend.numpoints'] = 1
    return plt

def parse_args(requires_snapshot=True, requires_observations=True):

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', help='SHArk configuration file')
    parser.add_argument('-m', '--model', help='Model name')
    parser.add_argument('-s', '--simu', help='Simulation name')
    parser.add_argument('-S', '--shark-dir', help='SHArk base output directory')
    parser.add_argument('-v', '--subvolumes', help='Comma- and dash-separated list of subvolumes to process', default='0')
    parser.add_argument('-o', '--output-dir', help='Output directory for plots. Defaults to <shark-dir>/Plots/<simu>/<model>')

    if requires_observations:
        parser.add_argument('-O ', '--obs-dir', help='Observations directory')

    if requires_snapshot:
        parser.add_argument('snapshot', help='Snapshot output to process', type=int)

    opts = parser.parse_args()
    if not opts.config and (not opts.model or not opts.simu or not opts.shark_dir):
        parser.error('Either -c or -m/-s/-S must be given')

    if requires_snapshot and opts.snapshot is None:
        parser.error('snapshot option is required')
    if requires_observations and opts.obs_dir is None:
        parser.error('-O is required')

    if opts.config:
        cparser = configparser.ConfigParser()
        cparser.read(opts.config)
        model = cparser.get('execution', 'name_model')
        simu = cparser.get('simulation', 'sim_name')
        shark_dir = cparser.get('execution', 'output_directory')
        print("Parsed configuration file %s" % (opts.config,))
    else:
        model = opts.model
        simu  = opts.simu
        shark_dir = opts.shark_dir
    model_dir = os.path.join(shark_dir, simu, model)

    output_dir = opts.output_dir
    if not output_dir:
        output_dir = os.path.join(shark_dir, 'Plots', simu, model)
    print("Creating plots under %s" % (output_dir,))

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    subvolumes = []
    for r in filter(None, opts.subvolumes.split(',')):
        if '-' in r:
            x = [int(x) for x in r.split('-')]
            subvolumes.extend(list(range(x[0], x[1])))
            subvolumes.append(x[1])
        else:
            subvolumes.append(int(r))

    ret = [model_dir, output_dir, subvolumes]
    if requires_observations:
        ret.append(opts.obs_dir)
    if requires_snapshot:
        ret.append(opts.snapshot)
    return ret

def load_observation(obsdir, fname, cols):
    fname = os.path.join(obsdir, fname)
    print("Loading observations from %s" % fname)
    return np.loadtxt(fname, usecols=cols, unpack=True)

def prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1)):

    from matplotlib.ticker import MultipleLocator

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    if xtit:
        ax.set_xlabel(xtit,fontsize=13)
    if ytit:
        ax.set_ylabel(ytit,fontsize=13)
    ax.xaxis.set_minor_locator(MultipleLocator(locators[0]))
    ax.xaxis.set_major_locator(MultipleLocator(locators[1]))
    ax.yaxis.set_minor_locator(MultipleLocator(locators[2]))
    if len(locators) > 3:
        ax.yaxis.set_major_locator(MultipleLocator(locators[3]))
    ax.tick_params(labelsize=12)

def prepare_legend(ax, colors, loc=None, **legend_kwargs):

    if(loc == None):
	loc = 3
    leg = ax.legend(loc=loc, prop={'size': 12}, **legend_kwargs)
    for color,text in zip(colors,leg.get_texts()):
        text.set_color(color)
    leg.draw_frame(False)

def errorbars(ax, x, y, yerrdn, yerrup, color, marker,
              err_absolute=True, condition=None, **kwargs):

    if condition is not None:
        ind = np.where(condition)
        x = x[ind]
        y = y[ind]
        yerrup = yerrup[ind]
        yerrdn = yerrdn[ind]

    if err_absolute:
        yerr = [y - yerrdn, yerrup - y]
    else:
        yerr = [yerrdn, yerrup]

    ax.errorbar(x, y, yerr=yerr, ls='None', mfc='None', ecolor=color, mec=color,
                marker=marker, **kwargs)

def savefig(output_dir, fig, plotname):
    plotfile = os.path.join(output_dir, plotname)
    print('Saving plot to %s' % plotfile)
    fig.savefig(plotfile, dvi=300, pad_inches=0)

def read_data(model_dir, snapshot, fields, subvolumes, include_h0_volh=True):
    """Read the galaxies.hdf5 file for the given model/snapshot/subvolume"""

    data = collections.OrderedDict()
    for idx, subv in enumerate(subvolumes):

        fname = os.path.join(model_dir, str(snapshot), str(subv), 'galaxies.hdf5')
        print('Reading data from %s' % fname)
        with h5py.File(fname, 'r') as f:
            if idx == 0 and include_h0_volh:
                data['h0'] = f['Cosmology/h'].value
                data['vol'] = f['runInfo/EffectiveVolume'].value * len(subvolumes)

            for gname, dsnames in fields.items():
                group = f[gname]
                for dsname in dsnames:
                    full_name = '%s/%s' % (gname, dsname)
                    l = data.get(full_name, None)
                    if l is None:
                        l = group[dsname].value
                    else:
                        l = np.concatenate([l, group[dsname].value])
                    data[full_name] = l

    return list(data.values())
