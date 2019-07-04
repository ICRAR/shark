#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2019
# Copyright by UWA (in the framework of the ICRAR)
#
# Originally contributed by Mawson Sammons
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
"""
This script produces the diagnostic plots used to help visualise the PSO process
One is a 3D representation of the swarm movement in the parameter space the other
is the Log Liklihood over iteration for each particle
"""


import argparse
import os
import re

import matplotlib
import matplotlib.animation as anim
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import analysis

try:
    import pandas
    import seaborn
except:
    seaborn = pandas = None


pos_re = re.compile('track_[0-9]+_pos.npy')
fx_re = re.compile('track_[0-9]+_fx.npy')

def load_space_and_particles(tracks_dir, space_file):
    """Loads the PSO pos/fx information stored within a directory, plus the
    original search space file"""

    all_fnames = list(os.listdir(tracks_dir))
    pos_fnames = list(filter(lambda x: pos_re.match(x), all_fnames))
    fx_fnames = list(filter(lambda x: fx_re.match(x), all_fnames))

    if len(pos_fnames) != len(fx_fnames):
        print("Different number of pos/fx files, using common files only")
        l = min(len(pos_fnames), len(fx_fnames))
        del pos_fnames[l:]
        del fx_fnames[l:]

    print("Loading %d pos/fx files" % len(pos_fnames))

    # Read files in filename order and populate pos/fx np arrays
    pos_fnames.sort()
    fx_fnames.sort()
    pos = []
    fx = []
    for pos_fname, fx_fname in zip(pos_fnames, fx_fnames):
        pos.append(np.load(os.path.join(tracks_dir, pos_fname)))
        fx.append(np.load(os.path.join(tracks_dir, fx_fname)))

    # after this fx dims are (S, L), pos dims are (S, D, L)
    pos, fx = np.asarray(pos), np.asarray(fx)
    pos = np.moveaxis(pos, 0, -1)
    fx = np.moveaxis(fx, 0, -1)

    space = analysis.load_space(space_file)
    if space.shape[0] != pos.shape[1]:
        raise ValueError("Particles have different dimensionality than space")

    return space, pos, fx

def plot_performance(fx, fig=None):
    """Creates a performance plot for all PSO particles across iterations"""

    ind = np.where(fx < 1e20)
    maxrange = max(fx[ind])
    minrange = min(fx[ind])
    ind = np.where(fx <= minrange)
    print ("particle with the smallest likelihood (particle number, iteration) %r" %(ind,))
    fig = fig or plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylim(minrange, maxrange)
    ax.set_ylabel('Likelihood evaluation')
    ax.set_xlabel('Iterations')
    for i, _fx in enumerate(fx):
        ax.plot(np.arange(_fx.shape[0]), _fx, label='Particle %d' % i)
    return fig

def plot_3d_space_evolution(space, pos, fx, fig=None):
    """Creates a 3D plot showing the evolution of PSO particles in a 3D search space"""

    cm = plt.get_cmap('plasma')
    cNorm = matplotlib.colors.Normalize(vmin=np.amin(fx), vmax=np.amax(fx))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

    fig = fig or plt.figure()
    ax = Axes3D(fig)
    for _pos, _fx in zip(pos, fx):
        ax.scatter(_pos[0,:], _pos[1,:], _pos[2,:], c=scalarMap.to_rgba(_fx[:]), s=(np.flip(np.arange(20), 0)/10)**9)

    scalarMap.set_array(fx)
    cbar=fig.colorbar(scalarMap)
    cbar.set_label('Log Likelihood', fontsize=16)

    axis_labels = space['plot_label']
    lb = space['lb']
    ub = space['ub']
    ax.set_xlabel(axis_labels[0], fontsize=16)
    ax.set_ylabel(axis_labels[1], fontsize=16)
    ax.set_zlabel(axis_labels[2], fontsize=16)
    ax.set_xlim(lb[0], ub[0])
    ax.set_ylim(lb[1], ub[1])
    ax.set_zlim(lb[2], ub[2])

    return fig

def plot_3d_space_animation(space, pos, fx, fig=None):
    """Creates a 3D plot showing the evolution of PSO particles in a 3D search space"""

    fig = fig or plt.figure()
    fig.set_tight_layout(True)

    # Color mapping setup
    cm = plt.get_cmap('plasma')
    cNorm = matplotlib.colors.Normalize(vmin=np.amin(fx), vmax=np.amax(fx))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    scalarMap.set_array(fx)

    # pos: S, D, L -> L, D, S / fx: S, L -> L, S
    pos = np.swapaxes(pos, 0, -1)
    fx = np.swapaxes(fx, 0, -1)

    axis_labels = space['plot_label']
    lb = space['lb']
    ub = space['ub']
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel(axis_labels[0])
    ax.set_ylabel(axis_labels[1])
    ax.set_zlabel(axis_labels[2])
    ax.set_xlim(lb[0], ub[0])
    ax.set_ylim(lb[1], ub[1])
    ax.set_zlim(lb[2], ub[2])
    scat = ax.scatter(pos[0][0], pos[0][1], pos[0][2],
                       c=scalarMap.to_rgba(fx[0]))
    title = ax.text2D(0.5, 0.95, 'Iteration 0', transform=fig.transFigure,
            ha='center', va='top')

    # Colorbar indicating color mapping
    cbar = fig.colorbar(scalarMap)
    cbar.set_label('Log Likelihood')
    colors = [scalarMap.to_rgba(x) for x in fx]

    def update(x):
        count, (pos, fx, color) = x
        scat._offsets3d = pos
        scat.set_color(color)
        title.set_text("Iteration %d" % (count + 1))
        return scat, title

    frames_data = list(enumerate(zip(pos, fx, colors)))
    animation = anim.FuncAnimation(fig, update, frames=frames_data, blit=False)
    return animation

def plot_pairplot(space, pos):
    """Produce a pairplot with seaborn"""
    S, D, L = pos.shape
    pos = np.swapaxes(pos, 0, 1)
    f = pandas.DataFrame(pos.reshape((D, S * L)).T, columns=space['plot_label'])
    return seaborn.pairplot(f)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-S', '--space-file',
        help='File with the search space specification, defaults to space.txt',
        default='space.txt')
    parser.add_argument('tracks_dir')
    opts = parser.parse_args()

    space, pos, fx = load_space_and_particles(opts.tracks_dir, opts.space_file)

    S, D, L = pos.shape
    print('Producing plots for S=%d, D=%d, L=%d' % (S, D, L))

    fig = plot_performance(fx)
    fig.savefig('performance.pdf')

    if D == 3:
        fig = plot_3d_space_evolution(space, pos, fx)
        fig.savefig('3DPSOC.pdf')
        animation = plot_3d_space_animation(space, pos, fx)
        animation.save('3d-particles.gif', writer='imagemagick', fps=8)

    if seaborn:
        fig = plot_pairplot(space, pos)
        fig.savefig('pairplot.pdf')

if __name__ == '__main__':
    main()
