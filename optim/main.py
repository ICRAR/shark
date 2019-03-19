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
import shutil
import subprocess
import sys
import tempfile
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
import numpy as np
import pso


logger = logging.getLogger('main')

if sys.version_info[0] == 3:
    b2s = lambda b: b.decode('ascii')
    raw_input = input
else:
    b2s = lambda b: b

def count_jobs(job_name):
    """Returns how many jobs with self.jobs_name are currently queued or running"""

    try:
        out, err, code = common.exec_command("squeue")
    except OSError:
        raise RuntimeError("Couldn't run squeue, is it installed?")

    if code:
        raise RuntimeError("squeue failed with code %d: stdout: %s, stderr: %s" % (code, out, err))

    lines_with_jobname = [l for l in out.splitlines() if job_name in l]
    return len(lines_with_jobname)

def _exec_shark(msg, cmdline):
    logger.info('%s with command line: %s', msg, subprocess.list2cmdline(cmdline))
    out, err, code = common.exec_command(cmdline)
    if code != 0:
        logger.error('Error while executing %s (exit code %d):\n' +
                     'stdout:\n%s\nstderr:\n%s', cmdline[0], code, b2s(out), b2s(err))
        raise RuntimeError('%s error' % cmdline[0])


def _to_shark_options(particle, space):
    """Given `particle` in `space` return an iterable with the corresponding
    shark options settings corresponding to that particle"""
    for value, name, is_log in zip(particle, space['name'], space['is_log']):
        if is_log:
            value = 10 ** value
        yield '%s=%s' % (name, value)


count = 0
def run_shark_hpc(particles, *args):
    """
    - Handler function for running PSO on Shark on a SLURM based cluster.
    - Swarm size and number of iterations need to be set within the script for now
    - Function needs the relative path to a Shark config file under the -c option
    - For now the subprocess call within must be altered if you are changing shark submit options
    - To find appropriate memory allocations peruse the initial output lines of each particle
    """

    global count

    opts, space, subvols, statTest = args

    # Prepare the file that will be used by the shark submission scripts
    # to determine which values shark will be run for. We put a final \n so the
    # final line gets properly counted by wc (used by shark-submit)
    shark_options = [
        ' '.join(['-o "%s"' % option for option in _to_shark_options(particle, space)])
        for particle in particles
    ]
    positions_fname = tempfile.mktemp('particle_positions.txt')
    logger.info('Creating particle positions file at %s', positions_fname)
    with open(positions_fname, 'wt') as f:
        f.write('\n'.join(shark_options) + '\n')

    # Submit the execution of multiple shark instances, one for each particle
    job_name = 'PSOSMF_%d' % count
    shark_output_base = os.path.join(opts.outdir, job_name)
    cmdline = ['./shark-submit', '-S', opts.shark_binary, '-w', opts.walltime,
               '-n', job_name, '-O', shark_output_base, '-E', positions_fname,
               '-V', ' '.join(map(str, subvols))]
    if opts.account:
        cmdline += ['-a', opts.account]
    if opts.queue:
        cmdline += ['-Q', opts.queue]
    if opts.nodes:
        cmdline += ['-N', str(opts.nodes)]
    else:
        cmdline += ['-m', opts.memory, '-c', str(opts.cpus)]
    cmdline.append(opts.config)
    _exec_shark('Queueing PSO particles', cmdline)

    # Actually wait for the jobs to finish...
    while count_jobs(job_name) > 0:
        time.sleep(10)

    ss = len(particles)
    fx = np.zeros([ss, 3])
    for i in range(ss):
        _, simu, model, _ = common.read_configuration(opts.config)
        particle_outdir = os.path.join(shark_output_base, str(i))
        modeldir = common.get_shark_output_dir(particle_outdir, simu, model)
        for j, constraint in enumerate(opts.constraints):
            y_obs, y_mod, err = constraint.get_data(modeldir, subvols)
            fx[i, j] = statTest(y_obs, y_mod, err)
        if not opts.keep:
            shutil.rmtree(particle_outdir)

    fx = np.sum(fx, 1)
    logger.info('Particles %r evaluated to %r', particles, fx)

    # this global count just tracks the number of iterations so they can be saved to different files
    count += 1

    return fx


def run_shark(particle, *args):

    opts, space, subvols, statTest = args

    pid = multiprocessing.current_process().pid
    shark_output_base = os.path.join(opts.outdir, 'output_%d' % pid)
    _, simu, model, _ = common.read_configuration(opts.config)
    modeldir = common.get_shark_output_dir(shark_output_base, simu, model)

    cmdline = [opts.shark_binary, opts.config,
               '-o', 'execution.output_directory=%s' % shark_output_base,
               '-o', 'execution.simulation_batches=%s' % ' '.join(map(str, subvols))]
    for option in _to_shark_options(particle, space):
        cmdline += ['-o', option]
    _exec_shark('Executing shark instance', cmdline)

    total = 0
    for constraint in opts.constraints:
        y_obs, y_mod, err = constraint.get_data(modeldir, subvols)
        total += statTest(y_obs, y_mod, err)

    logger.info('Particle %r evaluated to %f', particle, total)

    if not opts.keep:
        shutil.rmtree(shark_output_base)

    return total


def setup_logging(outdir):
    # Setup the logging
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
        f = run_shark_hpc
    else:
        n_cpus = multiprocessing.cpu_count()
        procs = min(n_cpus, ss)
        f = run_shark

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
        logger.info('%10s [%.1f - %.1f]' % (c.__class__.__name__, c.domain[0], c.domain[1]))
    logger.info('HPC mode: %d', opts.hpc_mode)
    if opts.hpc_mode:
        logger.info('    Account used to submit: %s', opts.account if opts.account else '')
        logger.info('    Queue to submit: %s', opts.queue if opts.queue else '')
        logger.info('    Walltime per submission: %s', opts.walltime)
        logger.info('    CPUs per instance: %d', opts.cpus)
        logger.info('    Memory per instance: %s', opts.memory)
        logger.info('    Nodes to use: %s', opts.nodes)

    while True:
        answer = raw_input('\nAre these parameters correct? (Yes/no): ')
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