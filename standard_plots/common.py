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

def _select_closest(i1, i2, z, redshifts):
    i1 = min(i1, len(redshifts) - 1)
    i2 = min(i2, len(redshifts) - 1)
    _z1 = redshifts[i1]
    _z2 = redshifts[i2]
    if abs(_z1 - z) < abs(_z2 - z):
        return i1
    return i2
_select_closest = np.vectorize(_select_closest, excluded={3})

class _redshift_table(object):
    """A class holding a redshift table. It calculates the snapshot for a given z"""

    def __init__(self, fname):

        table = np.loadtxt(fname, dtype={'formats': ('i4', 'f8'), 'names': ('snapshot', 'z')})
        snapshots = table['snapshot']
        z = table['z']

        # ensure snapshots are increasing, and redshifts decreasing
        if not np.all(snapshots[:-1] < snapshots[1:]):
            raise ValueError('Snapshots are not always increasing')
        if not np.all(z[:-1] > z[1:]):
            raise ValueError('Redshifts are not always decreasing')

        # Revert both; we want redshifts to be in ascending order
        self.z = z[::-1]
        self.snapshots = snapshots[::-1]

    def __getitem__(self, z):
        """Get the corresponding snapshots for the given redshift"""

        # We find the upper bound index using searchsorted, but then verify
        # which end of the range [z[i-1], z[i]] is closest to the requested z
        # value and use the corresponding index
        if not np.isscalar(z) and len(z) == 0:
            return []
        idx = np.searchsorted(self.z, z)
        prev_idx = np.maximum(idx - 1, 0)
        idx = _select_closest(idx, prev_idx, z, self.z)
        return self.snapshots[idx]

def load_matplotlib():
    import matplotlib.pyplot as plt
    plt.rcParams['legend.numpoints'] = 1
    return plt

def get_output_dir(shark_dir, simu, model):
    return os.path.join(shark_dir, 'Plots', simu, model)

def read_configuration(config):
    cparser = configparser.ConfigParser()
    cparser.read(config)
    shark_dir = cparser.get('execution', 'output_directory')
    model = cparser.get('execution', 'name_model')
    simu = cparser.get('simulation', 'sim_name')
    redshift_file = cparser.get('simulation', 'redshift_file')
    return shark_dir, simu, model, redshift_file

def parse_args(requires_observations=True):

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', help='SHArk configuration file')
    parser.add_argument('-m', '--model', help='Model name')
    parser.add_argument('-s', '--simu', help='Simulation name')
    parser.add_argument('-S', '--shark-dir', help='SHArk base output directory')
    parser.add_argument('-z', '--redshift-file', help='Redshift table file for the simulation')
    parser.add_argument('-v', '--subvolumes', help='Comma- and dash-separated list of subvolumes to process', default='0')
    parser.add_argument('-o', '--output-dir', help='Output directory for plots. Defaults to <shark-dir>/Plots/<simu>/<model>')

    opts = parser.parse_args()

    # If a configuration file is not passed, then we need to have all the
    # individual options given to us
    if not opts.config and (not opts.model or not opts.simu or not opts.shark_dir or not opts.redshift_file):
        parser.error('Either -c or -m/-s/-S/-z must be given')

    opts.obs_dir = os.path.normpath(os.path.abspath(os.path.join(__file__, '..', '..', 'data')))
    if requires_observations and opts.obs_dir is None:
        parser.error('-O is required')

    # The following allows using opts.config to set all these preferences,
    # but also allows users to override any of them with the individual values
    # given via the command-line.
    if opts.config:
        shark_dir, simu, model, redshift_file = read_configuration(opts.config)
        print("Parsed configuration file %s" % (opts.config,))
    if opts.model:
        model = opts.model
    if opts.simu:
        simu = opts.simu
    if opts.shark_dir:
        shark_dir = opts.shark_dir
    if opts.redshift_file:
        redshift_file = opts.redshift_file
    model_dir = os.path.join(shark_dir, simu, model)

    output_dir = opts.output_dir
    if not output_dir:
        output_dir = get_output_dir(shark_dir, simu, model)
    print("Creating plots under %s" % (output_dir,))

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    # Having the replace(',', ' ') allows us to separate subvolumes by command
    # and/or space
    subvolumes = set()
    for r in filter(None, opts.subvolumes.replace(',', ' ').split()):
        if '-' in r:
            x = [int(x) for x in r.split('-')]
            subvolumes.update(range(x[0], x[1] + 1))
        else:
            subvolumes.add(int(r))
    print("Considering the following subvolumes: %s" % ' '.join([str(x) for x in subvolumes]))

    ret = [model_dir, output_dir, _redshift_table(redshift_file), tuple(subvolumes)]
    if requires_observations:
        ret.append(opts.obs_dir)
    return ret

def load_observation(obsdir, fname, cols):
    fname = os.path.join(obsdir, fname)
    print("Loading observations from %s" % fname)
    return np.loadtxt(fname, usecols=cols, unpack=True)

def prepare_ax(ax, xmin, xmax, ymin, ymax, xtit, ytit, locators=(1, 1, 1, 1), fontsize=13):

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
    loc = 3 if loc is None else loc
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
                data['h0'] = f['cosmology/h'].value
                data['vol'] = f['run_info/effective_volume'].value * len(subvolumes)

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

# If called as a program, print information taken from a configuration file
# This simple functionality is used by shark-submit to easily find out where
# the plots have been produced, and save us the trouble to re-implement it
if __name__ == '__main__':
    action = sys.argv[1]
    shark_dir, simu, model, redshift_file = read_configuration(sys.argv[2])
    if action == 'output_dir':
        print(get_output_dir(shark_dir, simu, model))
    elif action == 'snapshots':
        z = list(map(float, ' '.join(sys.argv[3:]).split()))
        table = _redshift_table(redshift_file)
        print(' '.join(map(str, table[z])))
