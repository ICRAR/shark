#
#    ICRAR - International Centre for Radio Astronomy Research
#    (c) UWA - The University of Western Australia, 2019
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

import logging
import os
import shutil
import subprocess
import tempfile
import time

import numpy as np

import common
import constraints
import functools


logger = logging.getLogger(__name__)


class AbortedByUser(Exception):
    """Signals that the execution was aborted by the user"""
    pass


def job_is_alive(job_id):
    """Returns whether `job_id` is still "alive" (i.e., in the queue) or not"""
    try:
        out, err, code = common.exec_command(['squeue', '--noheader', '-j', job_id])
        out = common.b2s(out)
    except OSError:
        raise RuntimeError("Couldn't run squeue, is it installed?")
    if code:
        raise RuntimeError("squeue failed with code %d: stdout: %s, stderr: %s" % (code, out, err))
    return len([l for l in out.splitlines() if job_id in l]) > 0


def cancel_job(job_id):
    """Cancels `job_id` by removing it from the queue"""
    logger.info('Cancelling job %s', job_id)
    try:
        out, err, code = common.exec_command(['scancel', job_id])
    except OSError:
        raise RuntimeError("Couldn't run scancel, is it installed?")
    if code:
        raise RuntimeError("scancel failed with code %d: stdout: %s, stderr: %s" % (code, out, err))

    tries = 0
    max_tries = 20
    logger.debug('Checking that job %s has been successfully cancelled', job_id)
    while job_is_alive(job_id) and tries < max_tries:
        time.sleep(1)
        tries += 1
    if tries == max_tries:
        logger.warning('Job ID %s is still alive, you will need to cancel it manually', job_id)
    else:
        logger.info('Job %s successfully cancelled', job_id)


def _exec_shark(msg, cmdline):
    logger.info('%s with command line: %s', msg, subprocess.list2cmdline(cmdline))
    out, err, code = common.exec_command(cmdline)
    if code != 0:
        logger.error('Error while executing %s (exit code %d):\n' +
                     'stdout:\n%s\nstderr:\n%s', cmdline[0], code,
                     common.b2s(out), common.b2s(err))
        raise RuntimeError('%s error' % cmdline[0])
    return common.b2s(out)


def _to_shark_options(particle, space):
    """Given `particle` in `space` return an iterable with the corresponding
    shark options settings corresponding to that particle"""
    for value, name, is_log in zip(particle, space['name'], space['is_log']):
        if is_log:
            value = 10 ** value
        yield '%s=%s' % (name, value)


count = 0


def _prepare_run(opts, particles):
    iteration_name = 'PSO_%d' % count
    shark_output_base = os.path.join(opts.outdir, iteration_name)
    os.makedirs(shark_output_base)
    np.save(os.path.join(shark_output_base, 'particles.npy'), particles)
    return iteration_name, shark_output_base


def _read_particle_results(particle_index, shark_output_base, opts, subvols, statTest):
    _, simu, model, _ = common.read_configuration(opts.config)
    particle_outdir = os.path.join(shark_output_base, str(particle_index))
    modeldir = common.get_shark_output_dir(particle_outdir, simu, model)
    _y_obs, y_mod, y_err, result = constraints.evaluate(opts.constraints, statTest, modeldir, subvols)
    if not opts.keep:
        shutil.rmtree(particle_outdir, ignore_errors=True)
    return _y_obs, y_mod, y_err, result

def _read_results(particles, mp_pool, shark_output_base, opts, subvols, statTest):

    logger.info("=================== Evaluating constraints")
    y_obs, y_mod, y_err, results = zip(*mp_pool.map(
        functools.partial(
            _read_particle_results,
            shark_output_base=shark_output_base,
            opts=opts,
            subvols=subvols,
            statTest=statTest
        ),
        range(len(particles))
    ))

    constraints.log_results(opts.constraints, results)
    try:
        np.save(os.path.join(shark_output_base, 'modelvals.npy'), y_mod)
        np.save(os.path.join(shark_output_base, 'modelerrorvals.npy'), y_err)
    except ValueError:
        logger.error("Couldn't save model values/errors because they don't all have the same shape. TO BE IMPLEMENTED IN THE FUTURE.")

    results = np.sum(results, axis=1)

    # this global count just tracks the number of iterations so they can be
    # saved to different files
    global count
    count += 1

    return results


def _run_shark_hpc(particles, shark_output_base, job_name, opts, space, subvols):
    """
    - Handler function for running PSO on Shark on a SLURM based cluster.
    - Swarm size and number of iterations need to be set within the script for now
    - Function needs the relative path to a Shark config file under the -c option
    - For now the subprocess call within must be altered if you are changing shark submit options
    - To find appropriate memory allocations peruse the initial output lines of each particle
    """

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
    out = _exec_shark('Queueing PSO particles', cmdline)

    # In the output of shark-submit there will be the Job ID we need to monitor
    job_id = None
    _job_id_msg = 'shark-submit: New job submitted with ID '
    for line in out.splitlines():
        if line.startswith(_job_id_msg):
            job_id = line.split(_job_id_msg)[-1]
            break
    if job_id is None:
        raise RuntimeError("Couldn't get the ID of the new submitted job, cannot continue")

    # Actually wait for the job to finish...
    try:
        while job_is_alive(job_id):
            time.sleep(10)
    except KeyboardInterrupt:
        cancel_job(job_id)
        raise AbortedByUser


def _run_single_local_shark(indexed_particle, shark_output_base, opts, space, subvols):
    index, particle = indexed_particle
    shark_output_dir = os.path.join(shark_output_base, str(index))
    cmdline = [opts.shark_binary, opts.config,
               '-o', 'execution.output_directory=%s' % shark_output_dir,
               '-o', 'execution.simulation_batches=%s' % ' '.join(map(str, subvols))]
    for option in _to_shark_options(particle, space):
        cmdline += ['-o', option]
    try:
        _exec_shark('Executing shark instance', cmdline)
    except KeyboardInterrupt:
        raise AbortedByUser


def _run_shark_local(particles, mp_pool, shark_output_base, opts, space, subvols):
    mp_pool.map(
        functools.partial(
            _run_single_local_shark,
            shark_output_base=shark_output_base,
            opts=opts,
            space=space,
            subvols=subvols
        ),
        enumerate(particles)
    )


def run_shark(particles, mp_pool, opts, space, subvols, statTest):
    """
    - Handler function for running PSO on Shark on a SLURM based cluster.
    - Swarm size and number of iterations need to be set within the script for now
    - Function needs the relative path to a Shark config file under the -c option
    - For now the subprocess call within must be altered if you are changing shark submit options
    - To find appropriate memory allocations peruse the initial output lines of each particle
    """
    shark_output_base, instance_name = _prepare_run(opts, particles)
    if opts.hpc_mode:
        _run_shark_hpc(particles, shark_output_base, instance_name, opts, space, subvols)
    else:
        _run_shark_local(particles, mp_pool, shark_output_base, opts, space, subvols)
    return _read_results(particles, mp_pool, shark_output_base, opts, subvols, statTest)