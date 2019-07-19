Optimization
============

|s| ships with an optimization package
to allow users easily explore
different |s| parameter sets
in different scenarios.


.. _optim.search_space:

Specifying the search space
---------------------------

The parameter search space to be used by the optimization routine
is specified in a text file,
where each row fully describes a parameter.

Each row is a comma-separated list of values looking like this::

 parameter-name, plot-label, is_log, lb, up

where:

 * ``parameter-name`` is the name of the |s| option
   as given in the :ref:`command line <config.cmdline>`
   (e.g., ``reincorporation.tau_reinc``).
 * ``plot-label`` is the label that will appear
   on any plots produced to diagnose the execution.
 * ``is_log`` should be ``1`` if the parameter
   should be optimized for in the log space,
   or ``0`` if it should be optimized linearly.
 * ``lb`` and ``ub`` define the parameter's search space boundaries.
   If ``is_log`` is ``1`` these values should already be logged.


As an example::

 reincorporation.tau_reinc, $\log_{10}(\tau_{reinc})$, 1, 0, 1.4771212547196624
 reincorporation.mhalo_norm, $\log_{10}(M_{norm})$, 1, 9, 12
 reincorporation.halo_mass_power, $\gamma$, 0, -3, 0


.. _optim.running:

Running
-------

To run the optimization routines for |s|
you need to enter to the ``optim`` subdirectory
of the |s| git repository
and run the ``main.py`` script.

Run ``main.py -h`` to get a comprehensive list
of all options that can be passed down to the script,
including options specific to
:ref:`PSO <optim.methods.pso>`,
:ref:`HPC environments <optim.hpc>`
and more.


.. _optim.methods:

Optimization methods
--------------------

Currently the only supported optimization method
is the `Particle Swarm Optimization
<https://en.wikipedia.org/wiki/Particle_swarm_optimization>`_ (PSO),
but with time we plan to add more optimization methods.


.. _optim.methods.pso:

PSO
^^^

When running a PSO optimization for shark,
users can specify a number of options:

 * ``-s SWARM_SIZE`` indicates the swarm size
   (i.e., the number of particles to use).
   It defaults to ``10 + sqrt(D) * 2``,
   where ``D`` is the number of dimensions of the problem,
   which corresponds to the number of parameters being fitted.
 * ``-m MAX_ITERATIONS`` is the number of maximum iterations
   that PSO should run for before giving up.
   Otherwise PSO will automatically stop
   when the particles start converging within certain limits
   (``1e-8`` in particle step differences or objective function changes).


.. _optim.eval_funcs:

Evaluation functions
--------------------

When comparing model data against observational data,
different functions can be used
to evaluate how well they compare to each other.
We currently support two evaluation functions,
which are specified using the ``-t`` flag:

 * ``chi2``: A :math:`\chi^2` distribution
 * ``student-t``: A log Student-T distribution

In both cases the evaluation function is applied
to individual data points pairs (observations v/s model data),
and a final sum is done to get the final result.
Thus, constraints with more data points
have naturally more relevance in the final result.


.. _optim.constraints:

Constraints
-----------

A flexible number of constrains are supported
to evaluate how well a |s| parameter set behaves.
Constraints are specified on the command line
with the ``-x`` switch (see ``main.py -h`` for details).
Each constraint specification follows this pattern::

 <name>[(<min>-<max>)][*<weight>]

Here sections within ``[]`` are optional,
meaning that only ``<name>`` is required.
``<name>`` is the name of the constraint (see below),
``<min>`` and ``<max>`` specify the domain to consider
during the evaluation of the constraint,
and ``<weight>`` is the relative weight of the constraint
when evaluating in a multi-constraint scenario.
Each constraint has a hard-coded domain that is used as default,
in the user doesn't specify one.
``<weight>`` defaults to 1,
meaning that the results of the evaluation function for all constraints
(see :ref:`above <optim.eval_funcs>` for details on this)
are weighted equally.

The following constraints are currently supported:

 * ``HIMF``: it evaluates the HI mass function at ``z=0``.
   Its default domain is ``(7, 12)``.
 * ``SMF_z0``: it evaluates the stellar mass function at ``z=0``.
   Its default domain is ``(8, 13)``.
 * ``SMF_z1``: like ``SMF_z0`` but at ``z=1``.

If you want to add new constraints
please refer to the
:doc:`in depth documentation <constraints_impl>`
about the subject.


.. _optim.hpc:

HPC support
-----------

HPC support for running |s| optimizations
comes out of the box.
In particular,
the :ref:`parallel parameter set evaluation <hpc.param_sets>`
supported by the |ss| script is used
to execute as many |s| instances in parallel as possible
in the cluster.

To turn on HPC support use the ``-H`` option
when calling ``main.py``.
Use ``main.py -h`` to get a full list of parameters.
Remember that some of these can already be defined
via :ref:`environment variables <hpc.envvars>`,
easing the usage of the system.


.. _optim.diagnostics:

Diagnostics
-----------

After running,
the optimization routines will generate a series of files
under a ``tracks`` folder.
These can be visually analyzed by running the ``diagnostics.py`` script
pointing to the ``tracks`` folder.
