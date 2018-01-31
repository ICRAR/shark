# SHArk

SHArk (Semi-analytic Halo Ark)
is a new, flexible semi-analytic model of galaxy formation.

## Requirements

SHArk is written in C++11,
and therefore an up-to-date compiler supporting this standard is required.
Just as a reference, `gcc` supports C++11 since version 4.8.1
(although previous versions supported it partially),
and `clang` supports the standard since version 3.3.

SHArk depends on the following libraries:

 * [GSL](https://www.gnu.org/software/gsl/) >= 2.0
 * [HDF5](https://support.hdfgroup.org/HDF5/)
 * [Boost](http://www.boost.org/)

These libraries are usually available as packages
in most Linux distributions and MacOS package managers.
For example:

* In Debian/Ubuntu (as root): `apt-get install libgsl-dev libhdf5-dev libboost-dev`
* In Fedora/CentOS/RedHat (as root): `yum install gsl-devel hdf5-devel boost-devel`
* In MacOS + Homebrew: `brew install gsl hdf5`
* In MacOS + MacPorts: `port install gsl-devel && port install hdf5`

If compiling in a supercomputing facility,
it is likely that these libraries will be available as loadable modules.
Use `module avail xyz` to look for them (e.g., `module avail boost`),
then load them with `module load xyz` (e.g., `module load boost`).

## Compiling

SHArk uses [CMake](https://cmake.org/) as its build tool.
`cmake` is used to perform system-level checks,
like looking for libraries and setting up the rules for the build,
and then generates the actual build scripts
in one of the supported build systems.
Among other things, `cmake` supports out-of-tree builds
(useful to keep more than one build with different settings,
and to avoid cluttering the original source code directories)
and several build system, like `make` and `ninja` files.

Standing on the root your repository,
you can run `cmake`  to produce Makefiles
and then compile SHArk with these steps:

```
mkdir build
cd build
cmake .. # By default will generate Makefiles
make all
```

This will produce all executables and libraries
contained in this project
