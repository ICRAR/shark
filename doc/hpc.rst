HPC systems
===========

|s| can run not only in desktop computers and laptops
but also in High Performance Computing (HPC) clusters.
Some special considerations need to be taken into account
when using |s| in these environments,
which are detailed in this section.

In this section we first briefly explain
some basic HPC concepts
for people who are not familiar with HPC systems
(feel free to skip otherwise),
and then describe how to
:ref:`build <hpc.building>` and :ref:`run <hpc.running>` |s|
in these environments.


Brief introduction to HPC
-------------------------

.. note::
 Feel free to skip this section if you know
 your way around HPC systems already

Modules
^^^^^^^

One of the main differences between personal computers and HPC systems
is that in these environments
requirements (sometimes including the compiler itself)
are usually installed in the form of *modules*.
Modules are software installations
that can be loaded and unloaded interactively,
making the installed software within them
available and unavailable at will.
For a full description of how modules work
see `the modules documentation <https://modules.readthedocs.io/en/stable/index.html>`_.

Please refer to your HPC centre's documentation
if you need more details on installed modules,
or if you are missing a required dependency.

Queues
^^^^^^

.. note::
 Feel free to skip this section if you know
 what queues are and how they operate

HPC centres usually do not allow
executing programs in nodes interactively.
Instead, all users need to submit a *job* to a *queue*.
A job is simply a single script
containing all the steps needed to run your work
in an automatic way.
When submitting jobs,
users need to specify how many resources it will need
(e.g., how many CPUs, computers or memory it requires).
Jobs are then eventually taken out of the queue
and executed in the cluster nodes.
*When* the job is executed depends on a few factors,
like how many resources your job is requesting,
the amount of currently available resources,
and your priority, among others.

Several queueing systems exist,
with `SLURM <https://slurm.schedmd.com/>`_
and `PBS/TORQUE <http://www.adaptivecomputing.com/products/torque/>`_
being the two most popular ones.
Although from a high-level perspective
they both offer similar functionality,
the particulars on how to use them
are quite different.
For example, submitting jobs to the queue
is done with different commands on each system,
and resource requirements are specified differently.


.. _hpc.building:

Building
--------

Building |s| in HPC systems follows the same procedure
used to build it in personal computers and laptops,
with the additional task
of having to load the necessary modules
containing :ref:`the requirements <building.reqs>` for building |s|.
To do this, use ``module avail xyz`` to look for them
(e.g., ``module avail boost``)
then load them with ``module load xyz``
(e.g., ``module load boost``).
If all goes well,
``cmake`` will be able to find the requirements.

If you are missing a module you have a few alternatives:

* Contact your HPC cluster support and ask them to install the missing module
* Build the software yourself (and optionally install it as a personal module)

.. _hpc.building.intelcc:

Intel compiler
^^^^^^^^^^^^^^

Depending on the version,
OpenMP support on the Intel compiler is a bit difficult to identify.
Until version 3.9 ``cmake`` was not able to identify OpenMP support
for newer Intel compilers,
and simply using a version of ``cmake`` >= 3.9 will solve the issue.
If you find yourself in this situation,
a big warning message will appear when running ``cmake``
to alert you and guide you in what to do.


.. _hpc.running:

Running with |ss|
-----------------

As stated in :ref:`running.scalability`,
|s| scales thanks to its input data
being independently separated in *sub-volumes*,
and thanks to its multithreading support via OpenMP.
When running on an HPC cluster
one should take advantage of these two aspects
to efficiently execute |s| across the available resources.

To this end |s| ships with a |ss| script
that can be used to easily submit jobs
that will run |s| over a list of sub-volumes
in an HPC cluster.
Internally, the script will spawn independent |s| instances
for each of the requested sub-volumes,
and will instruct each instance to use multiple CPUs.

|ss| not only eases this submission process;
it also abstracts away the differences between queueing systems,
making it easy for users to specify exactly what they want
without having to remember particular commands and formats.
It also creates the submissions in such a way
that all artifacts resulting from a submission
(i.e., log files, environment information, even plots)
end up in separate, well defined locations,
making it easy to inspect them
and shared them if needed.

.. note::
 Only SLURM is currently supported,
 but TORQUE will be supported,
 and more systems might come with time.

The |ss| script lives
under the ``hpc`` subdirectory of the |s| repository.
Its most basic usage looks like this::

 $> hpc/shark-submit -V 0 <config_file>

That will submit an execution of |s| only for sub-volume 0
using the ``config_file`` configuration file.


.. _hpc.param_sets:

Running with different parameter sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As explained above,
the default mode of execution of |ss|
parallelizes |s| executions by sub-volume,
using the same configuration.
However, a second mode is supported,
where users can specify different parameter sets
to be evaluated against the same inputs.
This is useful, for instance,
when one is sampling a parameter search space
to :doc:`optimize shark <optim>` against certain constraints,
or during other exploratory exercises.

This mode is triggered by using the ``-E file`` flag.
``file`` must contain all the command-line flags
that will be given to each |s| instance, one per row.
For example::

 -o "reincorporation.tau_reinc=9.789522070051014" -o "reincorporation.mhalo_norm=53167281575.647736" -o "reincorporation.halo_mass_power=-1.662049864221243"
 -o "reincorporation.tau_reinc=4.433571656151598" -o "reincorporation.mhalo_norm=344951728442.8235" -o "reincorporation.halo_mass_power=-2.197944428980997"
 -o "reincorporation.tau_reinc=8.744659838237162" -o "reincorporation.mhalo_norm=106081569566.5114" -o "reincorporation.halo_mass_power=-2.2146743876637798"
 -o "reincorporation.tau_reinc=5.568250109183069" -o "reincorporation.mhalo_norm=36854502778.199234" -o "reincorporation.halo_mass_power=-1.909464612909543"
 -o "reincorporation.tau_reinc=2.9521943986079364" -o "reincorporation.mhalo_norm=1916185243.0129645" -o "reincorporation.halo_mass_power=-2.797548035509215"

When this mode is used,
the ``-V`` flag that indicates subvolumes
applies equally to all |s| instances.
For example, if the user specifies ``-V 0-3 -E simple.txt``,
and ``simple.txt`` contains three lines,
then only three |s| instances will be spawned.



Options
^^^^^^^

|ss| supports many options,
which are roughly grouped into the following categories:

* *Queueing*: they include which queue to submit to,
  how many resources are needed (memory, CPUs and/or nodes),
  and more.
* *Plotting*: these control whether to produce
  the standard plots.
* *Shark*: these are |s|-specific options, like
  which particular |s| binary to use,
  and which sub-volumes to process and the configuration file to use.
* *Other*: modules to load, output directory to use,
  etc.

For a full help on all available options run::

 $> hpc/shark-submit -h

.. _hpc.envvars:

Environment variables
^^^^^^^^^^^^^^^^^^^^^

Some of the options of |ss|
will probably remain the same
across most (if not all) executions.
Because of these, a handful of environment variables
are inspected by |ss| and interpreted
as the default value for some of these options
(run ``shark-submit -h`` for a full list).
You can thus define these variables once
(e.g., in your ``~/.bash_rc`` or ``~/.bash_profile`` files)
to avoid having to repeat typing them each time.
