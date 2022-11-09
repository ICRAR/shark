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

import functools
import os
import multiprocessing

import common
import coldgas
import global_quantities
import hmf
import hothalo
import sizes
import smf
import smhm
import qso_lf

def main():

    model_dir, output_dir, redshift_table, subvolumes, obs_dir = common.parse_args()
    args_minimal = (model_dir, output_dir, redshift_table, subvolumes)
    args_with_obsdir = args_minimal + (obs_dir,)

    # Modules and which arguments they take
    args_and_mods = {
        args_minimal: [hothalo],
        args_with_obsdir: (coldgas, global_quantities, hmf, sizes, smf, smhm, qso_lf),
    }

    n_mods = functools.reduce(lambda x, y: x + y, [len(l) for l in args_and_mods.values()])
    n_procs = int(os.environ.get('SHARK_PLOT_PROCS', n_mods))
    print("Using %d processes to produce all plots" % n_procs)
    pool = multiprocessing.Pool(n_procs)

    # Go, go, go!
    futures = []
    for args, mods in args_and_mods.items():
        futures += [pool.apply_async(m.main, args) for m in mods]

    # Wait for all results to finish
    for f in futures:
        f.get()

if __name__ == '__main__':
    main()
