Running
=======

.. _running.cmdline:

Command line usage
------------------

.. note::
 These instructions apply for manual runs of |s|.
 If you want to run |s| in an HPC system
 please refer to :ref:`hpc.running`.

After successfully :doc:`building` |s|
you can directly run it from the command line::

 $> ./shark -h

This will print out detailed information
about how to run |s|.

|s| requires, and accepts, a number of command line options:

 * ``-h`` or ``-?`` show the help message and exits.
 * ``-V`` shows version information and exits.
 * ``-v <verbosity>`` specifies how verbose should |s| be.
   Values for ``verbosity`` range from 0 to 5,
   with 0 being mute and 5 being extremely verbose.
 * ``-t <threads>`` specifies how many OpenMP threads
   to use to run |s|.
   If 0 is given, then OpenMP will use its own default value.
   If |s| is compiled without OpenMP support
   this option is ignored.
 * ``-o <option>`` specifies additional configuration values
   to use. See :doc:`configuration/specifying` for details.

Any other argument is interpreted
as the name of a configuration file to load.
At least one configuration file is required for |s| to run.
See :doc:`configuration/specifying` for details on how configuration works.


Exit code
---------

Upon a successful run,
|s| returns with an exit code equals to ``0``,
or something different from ``0`` in case of any error.


.. _running.scalability:

Scalability
-----------

|s| scales using two approaches:

* Using **independent input data**, and
* Using **multiple CPUs** with OpenMP.

We describe both approached here,
and mention when they should be used.

Input data
----------

|s| input data (*volumes*) is usually divided
into separate files (*sub-volumes*).
These sub-volumes are self-sufficient,
meaning that they can be processed independently.

Based on this,
a |s| instance can be commanded to process
one or more sub-volumes.
This extremely simple but flexible scheme
allows for easy parallelisation
based on input data.
In other words,
multiple |s| instances can be independently executed
to process a big number of sub-volumes in parallel.
Note that this strategy
does not require communication between instances,
reducing both the complexity of the software and its dependencies
(i.e., it does not require MPI to work).

This is the basic strategy used by the |ss| script
when running under an :ref:`HPC environment <hpc.running>`.

OpenMP
------

During the main |s| evolution loop,
and for any snapshot,
the evolution of galaxies belonging to different merger trees
is independent from each other.
This is the place in the code
in which most of the time is spent,
and thus |s| parallelises the evolution of individual merger trees
so they take place in different threads.
Other parts of the code are parallelised as well.

|s| uses OpenMP to carry out this parallelisation.
OpenMP is widely supported by compilers nowadays,
but not universally (see :doc:`building` for details).
The number of threads to use is
specified on the command-line
(see :ref:`running.cmdline`, option ``-t``),
and can be set to either a fixed number,
or to the default value provided by the OpenMP library.

Using OpenMP will result in a speed up in most cases,
but only up to certain threshold
when using more CPUs will not necessarily improve
the runtime of |s|.
