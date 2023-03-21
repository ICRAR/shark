# Adapted from pyswarm verison 0.7 for optimising the shark SAM
#
# https://github.com/tisimst/pyswarm/tree/master/pyswarm   (original pyswarms code)

import logging
import os
import pickle
from functools import partial
import numpy as np


logger = logging.getLogger(__name__)

def _obj_wrapper(func, args, kwargs, x, iteration):
    return func(x, iteration, *args, **kwargs)

def _is_feasible_wrapper(func, x):
    return np.all(func(x)>=0)

def _cons_none_wrapper(x):
    return np.array([0])

def _cons_ieqcons_wrapper(ieqcons, args, kwargs, x):
    return np.array([y(x, *args, **kwargs) for y in ieqcons])

def _cons_f_ieqcons_wrapper(f_ieqcons, args, kwargs, x):
    return np.array(f_ieqcons(x, *args, **kwargs))

def _update_particle_pos_and_vel(pos_lb, pos_ub, omega, phip, phig, pos, vel, best_pos, best_global_pos):
    rp = np.random.uniform(size=pos.shape)
    rg = np.random.uniform(size=pos.shape)

    # Update the particles velocities
    vel = omega * vel + phip * rp * (best_pos - pos) + phig * rg * (best_global_pos - pos)

    # Update the particles' positions
    pos = pos + vel
    # Correct for bound violations
    maskl = pos < pos_lb
    masku = pos > pos_ub
    pos = pos * (~np.logical_or(maskl, masku)) + pos_lb * maskl + pos_ub * masku
    return pos, vel


def _dump_pso_state(state_filename, *pso_state):
    logger.info("Storing PSO state into %s", state_filename)
    with open(state_filename, 'wb') as f:
        pickle.dump(pso_state, f)

def _load_pso_state(state_filename):
    logger.info("Loading PSO state from %s", state_filename)
    with open(state_filename, 'rb') as f:
        return pickle.load(f)

def pso(func, lb, ub, ieqcons=[], f_ieqcons=None, args=(), kwargs={}, 
        swarmsize=100, omega=0.5, phip=0.5, phig=0.5, maxiter=100, 
        minstep=1e-8, minfunc=1e-8, processes=1,
        dumpfile_prefix=None, state_filename=None):
    """
    Perform a particle swarm optimization (PSO)
   
    Parameters
    ==========
    func : function
        The function to be minimized
    lb : array
        The lower bounds of the design variable(s)
    ub : array
        The upper bounds of the design variable(s)
   
    Optional
    ========
    ieqcons : list
        A list of functions of length n such that ieqcons[j](x,*args) >= 0.0 in 
        a successfully optimized problem (Default: [])
    f_ieqcons : function
        Returns a 1-D array in which each element must be greater or equal 
        to 0.0 in a successfully optimized problem. If f_ieqcons is specified, 
        ieqcons is ignored (Default: None)
    args : tuple
        Additional arguments passed to objective and constraint functions
        (Default: empty tuple)
    kwargs : dict
        Additional keyword arguments passed to objective and constraint 
        functions (Default: empty dict)
    swarmsize : int
        The number of particles in the swarm (Default: 100)
    omega : scalar
        Particle velocity scaling factor (Default: 0.5)
    phip : scalar
        Scaling factor to search away from the particle's best known position
        (Default: 0.5)
    phig : scalar
        Scaling factor to search away from the swarm's best known position
        (Default: 0.5)
    maxiter : int
        The maximum number of iterations for the swarm to search (Default: 100)
    minstep : scalar
        The minimum stepsize of swarm's best position before the search
        terminates (Default: 1e-8)
    minfunc : scalar
        The minimum change of swarm's best objective value before the search
        terminates (Default: 1e-8)
    processes : int
        The number of processes to use to evaluate objective function and 
        constraints. If processes = 0 then all particles are given to a single
        handling function to deal with them all at once (default: 1)
   
    Returns
    =======
    g : array
        The swarm's best known position (optimal design)
    f : scalar
        The objective value at ``g``
    p : array
        The best known position per particle
    pf: arrray
        The objective values at each position in p
   
    """
   
    assert len(lb)==len(ub), 'Lower- and upper-bounds must be the same length'
    assert hasattr(func, '__call__'), 'Invalid function handle'
    lb = np.array(lb)
    ub = np.array(ub)
    assert np.all(ub>lb), 'All upper-bound values must be greater than lower-bound values'

    vhigh = np.abs(ub - lb)
    vlow = -vhigh

    # Initialize objective function
    obj = partial(_obj_wrapper, func, args, kwargs)
    update_particle_pos_and_vel = partial(_update_particle_pos_and_vel, lb, ub, omega, phip, phig)

    # Initialize dumping function if required
    if dumpfile_prefix:
        def dump(i, x, fx):
            np.save(dumpfile_prefix % i + "_fx", fx)
            np.save(dumpfile_prefix % i + "_pos", x)
    else:
        dump = lambda *_: None

    # Check for constraint function(s) #########################################
    if f_ieqcons is None:
        if not len(ieqcons):
            logger.debug('No constraints given')
            cons = _cons_none_wrapper
        else:
            logger.debug('Converting ieqcons to a single constraint function')
            cons = partial(_cons_ieqcons_wrapper, ieqcons, args, kwargs)
    else:
        logger.debug('Single constraint function given in f_ieqcons')
        cons = partial(_cons_f_ieqcons_wrapper, f_ieqcons, args, kwargs)
    is_feasible = partial(_is_feasible_wrapper, cons)

    # Initialize the multiprocessing module if necessary
    if processes > 1:
        import multiprocessing
        mp_pool = multiprocessing.Pool(processes)

    # Initialize the particle swarm ############################################
    if os.path.isfile(state_filename):
        first_iteration, x, v, p, fp, g, fg = _load_pso_state(state_filename)
    else:
        S = swarmsize
        D = len(lb)  # the number of dimensions each particle has
        p = np.zeros((S, D))  # best particle positions
        fp = np.ones(S) * np.inf  # best particle function values
        g = []  # best swarm position
        fg = np.inf  # best swarm position starting value
        first_iteration = 0
        x = None
        v = None

    # Iterate until termination criterion met ##################################
    for iteration in range(first_iteration, maxiter):

        # Get particle positions and velocities
        if x is None:
            assert v is None
            x = lb + np.random.rand(S, D) * (ub - lb)
            v = vlow + np.random.rand(S, D) * (vhigh - vlow)
        else:
            x, v = update_particle_pos_and_vel(x, v, p, g)

        # Update objectives and constraints
        if processes > 1:
            fx = np.array(mp_pool.map(obj, x, iteration))
            fs = np.array(mp_pool.map(is_feasible, x))
        elif processes != 0:
            for i in range(S):
                fx[i] = obj(x[i, :], iteration)
                fs[i] = is_feasible(x[i, :])
        else:
            fx = obj(x, iteration)
            fs = is_feasible(x)
        dump(iteration, x, fx)

        # Store particle's best position (if constraints are satisfied)
        i_update = np.logical_and((fx < fp), fs)
        p[i_update, :] = x[i_update, :].copy()
        fp[i_update] = fx[i_update]

        # Update swarm's best position
        i_min = np.argmin(fp)
        if fp[i_min] < fg:
            logger.info('New best for swarm at iteration %d: %r %.3f', iteration, p[i_min, :], fp[i_min])

            p_min = p[i_min, :].copy()
            stop = False
            if iteration != 0:
                if np.abs(fg - fp[i_min]) <= minfunc:
                    logger.info('Stopping search: Swarm best objective change less than %.3f', minfunc)
                    stop = True
                elif np.sqrt(np.sum((g - p_min)**2)) <= minstep:
                    logger.info('Stopping search: Swarm best position change less than %.3f', minstep)
                    stop = True
            g = p_min
            fg = fp[i_min]
            if stop:
                break
        elif iteration == 0:
            # At the start, there may not be any feasible starting point, so just
            # give iteration a temporary "best" point since iteration's likely to change
            g = x[0, :].copy()

        logger.info('Best after iteration %d: %r %.3f', iteration, g, fg)
        _dump_pso_state(state_filename, iteration + 1, x, v, p, fp, g, fg)
    else:
        logger.info('Stopping search: maximum iterations reached --> %d', maxiter)

    if not is_feasible(g):
        logger.warning("However, the optimization couldn't find a feasible design. sorry n.n")

    return g, fg, p, fp
