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
import argparse
import logging
import math
import multiprocessing
import os
import sys
import time

def _abspath(p):
    return os.path.normpath(os.path.abspath(p))

# TODO: this is only for convenience and to avoid copying code over
#
# Once these modules and those under standard_scripts are put together under
# a proper package structure this will be gone
sys.path.insert(0, _abspath(os.path.join(__file__, '..', '..', 'standard_plots')))

import analysis
import common
import constraints
import execution
import pso


logger = logging.getLogger('main')

def setup_logging(outdir):
    log_fname = os.path.join(outdir, 'shark_pso.log')
    fmt = '%(asctime)-15s %(name)s#%(funcName)s:%(lineno)s %(message)s'
    fmt = logging.Formatter(fmt)
    fmt.converter = time.gmtime
    logging.root.setLevel(logging.INFO)
    h = logging.StreamHandler(stream=sys.stdout)
    h.setFormatter(fmt)
    logging.root.addHandler(h)
    h = logging.FileHandler(log_fname)
    h.setFormatter(fmt)
    logging.root.addHandler(h)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', help='Configuration file used as the basis for running shark', type=_abspath)
    parser.add_argument('-v', '--subvolumes', help='Comma- and dash-separated list of subvolumes to process', default='0')
    parser.add_argument('-b', '--shark-binary', help='The shark binary to use, defaults to either "shark" or "../build/shark"',
                        default=None, type=_abspath)
    parser.add_argument('-o', '--outdir', help='Auxiliary output directory, defaults to .', default=_abspath('.'),
                        type=_abspath)
    parser.add_argument('-k', '--keep', help='Keep temporary output files', action='store_true')

    pso_opts = parser.add_argument_group('PSO options')
    pso_opts.add_argument('-s', '--swarm-size', help='Size of the particle swarm. Defaults to 10 + sqrt(D) * 2 (D=number of dimensions)',
                          type=int, default=None)
    pso_opts.add_argument('-m', '--max-iterations', help='Maximum number of iterations to reach before giving up, defaults to 20',
                          default=10, type=int)
    pso_opts.add_argument('-S', '--space-file', help='File with the search space specification, defaults to space.txt',
                          default='space.txt', type=_abspath)
    pso_opts.add_argument('-t', '--stat-test', help='Stat function used to calculate the value of a particle, defaults to student-t',
                          default='student-t', choices=list(analysis.stat_tests.keys()))
    pso_opts.add_argument('-x', '--constraints', default='HIMF,SMF_z0,SMF_z1',
                          help=("Comma-separated list of constraints, any of HIMF, SMF_z0 or SMF_z1, defaults to 'HIMF,SMF_z0,SMF_z1'. "
                                "Can specify a domain range after the name (e.g., 'SMF_z0(8-11)')"))

    hpc_opts = parser.add_argument_group('HPC options')
    hpc_opts.add_argument('-H', '--hpc-mode', help='Enable HPC mode', action='store_true')
    hpc_opts.add_argument('-C', '--cpus', help='Number of CPUs per shark instance', default=1, type=int)
    hpc_opts.add_argument('-M', '--memory', help='Memory needed by each shark instance', default='1500m')
    hpc_opts.add_argument('-N', '--nodes', help='Number of nodes to use', default=None, type=int)
    hpc_opts.add_argument('-a', '--account', help='Submit jobs using this account', default=None)
    hpc_opts.add_argument('-q', '--queue', help='Submit jobs to this queue', default=None)
    hpc_opts.add_argument('-w', '--walltime', help='Walltime for each submission, defaults to 1:00:00', default='1:00:00')

    opts = parser.parse_args()

    if not opts.config:
        parser.error('-c option is mandatory but missing')

    if opts.shark_binary and not common.has_program(opts.shark_binary):
        parser.error("shark binary '%s' not found, specify a correct one via -b" % opts.shark_binary)
    elif not opts.shark_binary:
        for candidate in ['shark', '../build/shark']:
            if not common.has_program(candidate):
                continue
            opts.shark_binary = _abspath(candidate)
            break
        if not opts.shark_binary:
            parser.error("No shark binary found, specify one via -b")

    _, _, _, redshift_file = common.read_configuration(opts.config)
    redshift_table = common._redshift_table(redshift_file)
    subvols = common.parse_subvolumes(opts.subvolumes)

    setup_logging(opts.outdir)

    opts.constraints = constraints.parse(opts.constraints)
    for c in opts.constraints:
        c.redshift_table = redshift_table

    # Read search space specification, which is a comma-separated multiline file,
    # each line containing the following elements:
    #
    # param_name, plot_label, is_log, lower_bound, upper_bound
    space = analysis.load_space(opts.space_file)

    ss = opts.swarm_size
    if ss is None:
        ss = 10 + int(2 * math.sqrt(len(space)))

    args = (opts, space, subvols, analysis.stat_tests[opts.stat_test])

    if opts.hpc_mode:
        procs = 0
        f = execution.run_shark_hpc
    else:
        n_cpus = multiprocessing.cpu_count()
        procs = min(n_cpus, ss)
        f = execution.run_shark

    logger.info('-----------------------------------------------------')
    logger.info('Runtime information')
    logger.info('    shark binary: %s', opts.shark_binary)
    logger.info('    Base configuration file: %s', opts.config)
    logger.info('    Subvolumes to use: %r', subvols)
    logger.info('    Output directory: %s', opts.outdir)
    logger.info('    Keep temporary output files: %d', opts.keep)
    logger.info("PSO information:")
    logger.info('    Search space parameters: %s', ' '.join(space['name']))
    logger.info('    Swarm size: %d', ss)
    logger.info('    Maximum iterations: %d', opts.max_iterations)
    logger.info('    Lower bounds: %r', space['lb'])
    logger.info('    Upper bounds: %r', space['ub'])
    logger.info('    Test function: %s', opts.stat_test)
    logger.info('Constraints:')
    for c in opts.constraints:
        logger.info('    %s', c)
    logger.info('HPC mode: %d', opts.hpc_mode)
    if opts.hpc_mode:
        logger.info('    Account used to submit: %s', opts.account if opts.account else '')
        logger.info('    Queue to submit: %s', opts.queue if opts.queue else '')
        logger.info('    Walltime per submission: %s', opts.walltime)
        logger.info('    CPUs per instance: %d', opts.cpus)
        logger.info('    Memory per instance: %s', opts.memory)
        logger.info('    Nodes to use: %s', opts.nodes)

    while True:
        answer = common.raw_input('\nAre these parameters correct? (Yes/no): ')
        if answer:
            if answer.lower() in ('n', 'no'):
                logger.info('Not starting PSO, check your configuration and try again')
                return
            print("Please answer 'yes' or 'no'")
            continue
        break

    # Directory where we store the intermediate results
    tracksdir = os.path.join(opts.outdir, 'tracks')
    try:
        os.makedirs(tracksdir)
    except OSError:
        pass

    # Go, go, go!
    logger.info('Starting PSO now')
    tStart = time.time()
    if opts.hpc_mode:
        os.chdir('../hpc')
    xopt, fopt = pso.pso(f, space['lb'], space['ub'], args=args, swarmsize=ss,
                         maxiter=opts.max_iterations, processes=procs,
                         dumpfile_prefix=os.path.join(tracksdir, 'track_%03d'))
    tEnd = time.time()

    global count
    logger.info('Number of iterations = %d', count)
    logger.info('xopt = %r', xopt)
    logger.info('fopt = %r', fopt)
    logger.info('PSO finished in %.3f [s]', tEnd - tStart)

if __name__ == '__main__':
    main()