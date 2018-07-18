Building
========

Compiler
--------

|s| is written in C++11,
so in principle any compiler and library supporting this standard can be used.
As a reference, ``gcc`` supports C++11 since version 4.8.1
(but its C++ library support did not come to shape until gcc 5.0)
and ``clang`` supports the standard since version 3.3.

.. _building.reqs:

Requirements
------------

|s| depends on the following libraries:

* `GSL <https://www.gnu.org/software/gsl/>`_ >= 2.0
* `HDF5 <https://support.hdfgroup.org/HDF5/>`_ >= 1.8.0
* `Boost <http://www.boost.org/>`_ >= 1.54

Optional requirements are:

* An `OpenMP <http://www.openmp.org/>`_-enabled compiler
* `CxxTest <https://cxxtest.com/>`_, required only to compile the unit tests

These libraries are usually available as packages
in most Linux distributions and MacOS package managers.
For example:

* In Debian/Ubuntu: ``sudo apt-get install libgsl-dev libhdf5-dev libboost-all-dev``
* In Fedora/CentOS/RedHat (as root): ``yum install gsl-devel hdf5-devel boost-devel``
* In MacOS + Homebrew: ``brew install gsl hdf5 boost``
* In MacOS + MacPorts: ``port install gsl-devel && port install hdf5 && port install boost``

If compiling in a supercomputing facility,
it is likely that these libraries will be available as loadable modules.
See :doc:`hpc` for details.

Compiling
---------

|s| uses `CMake <https://cmake.org/>`_ as its build tool.
``cmake`` is used to perform system-level checks,
like looking for libraries and setting up the rules for the build,
and then generates the actual build scripts
in one of the supported build systems.
Among other things, ``cmake`` supports out-of-tree builds
(useful to keep more than one build with different settings,
and to avoid cluttering the original source code directories)
and several build system, like ``make`` and ``ninja`` files.

Standing on the root your repository,
you can run ``cmake``  to produce Makefiles
and then compile |s| with these steps::

 $> mkdir build
 $> cd build
 $> cmake ..
 $> make all

Other tips:

* Use, for example, ``make all -j4`` to compile using 4 parallel tasks.
* Use ``make VERBOSE=1`` to see exactly what the compilation is doing
* Use ``make clean`` to remove all compiled code

To make sure your build was successful, run::

 $> ./shark -V

That should output the version of |s| plus other information.

Compilation options
^^^^^^^^^^^^^^^^^^^

With ``cmake`` you can also specify additional compilation flags
and compilation options via ``-Dname=value``.
``cmake`` comes with a number of built-in flags
that can be used to specify different options.
We present a few useful ones here,
but for a more comprehensive list please see
`cmake's help <https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/Useful-Variables>`_ on the topic,
or run ``cmake --help-variable-list``:

* ``CMAKE_BUILD_TYPE``: describes the type of build you are doing.
  It must be one of ``Debug``, ``Release``, ``RelWithDebInfo`` or ``MinSizeRel``.
  |s| defaults to ``Release``.
  Depending on the build type, different compilation options will be used
  to compile the source file. For example, in most compilers
  ``Release`` will add the ``-O3`` compilation option,
  while ``RelWithDebInfo`` will add ``-O2 -g``.
* ``CMAKE_CXX_COMPILER``: the C++ compiler to use.
  It is usually automatically detected, but you can force
  it to be a particular one.
* ``CMAKE_CXX_COMPILER``: additional C++ flags used to compile.
  These are applied on top of those defined by the ``CMAKE_BUILD_TYPE``.

|s| defines its own ``cmake`` flags:

* ``SHARK_TEST``: if ``ON`` it enables the compilation of unit tests.
* ``SHARK_NO_OPENMP``: if ``ON`` it disables OpenMP support.

Examples
^^^^^^^^

* Compile in debug mode using the ``clang++`` compiler, without OpenMP
  support::

   $> cmake .. -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Debug -DSHARK_NO_OPENMP=ON

* Compile tests together with the code::

   $> cmake .. -DSHARK_TEST=ON

* Compile against the GSL installation found under /opt/gsl::

   $> cmake .. -DGSL_ROOT_DIR=/opt/gsl

* Generate the fastest possible code for your local machine/architecture::

   $> cmake .. -DCMAKE_CXX_FLAGS="-march=native"
