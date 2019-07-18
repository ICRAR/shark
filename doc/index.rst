shark
=====

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   building
   configuration/index
   running
   hpc
   output_files
   optim
   constraints_impl
   changelog


|s| is a new, flexible semi-analytic model of galaxy formation.

Thanks to its flexibility,
|s| allows for easy exploration of different physical processes.
Shark has been implemented
with several models for gas cooling,
active galactic nuclei,
stellar and photo-ionization feedback,
and star formation (SF).
The software can determine
the stellar mass function and stellar–halo mass relation at z=0–4;
cosmic evolution of the star formation rate density,
stellar mass, atomic and molecular hydrogen;
local gas scaling relations;
and structural galaxy properties.
It performs particularly well
for the mass–size relation for discs/bulges,
the gas–stellar mass and stellar mass–metallicity relations.

|s| is written in C++11 and has been parallelized with OpenMP.
It currently compiles with all major compilers
(gcc, clang, MSVC), but any C++11-enabled compiler should work.
It comes with a set of standard plotting scripts, HPC-related utilities
to ease its usage across as many platforms as possible,
and optimization routines for easy parameter exploration.

Citing
------

If you are using |s| for your projects,
please cite the following paper,
which is the first one describing shark in full::

 @article{doi:10.1093/mnras/sty2440,
     author = {Lagos, Claudia del P and Tobar, Rodrigo J and Robotham, Aaron S G and Obreschkow, Danail and Mitchell, Peter D and Power, Chris and Elahi, Pascal J},
     title = {Shark: introducing an open source, free, and flexible semi-analytic model of galaxy formation},
     journal = {Monthly Notices of the Royal Astronomical Society},
     volume = {481},
     number = {3},
     pages = {3573-3603},
     year = {2018},
     doi = {10.1093/mnras/sty2440},
     URL = {http://dx.doi.org/10.1093/mnras/sty2440},
     eprint = {/oup/backfile/content_public/journal/mnras/481/3/10.1093_mnras_sty2440/1/sty2440.pdf}
 }

An online entry can also be found
at `NASA's ADS service <https://ui.adsabs.harvard.edu/?#abs/2018MNRAS.481.3573L/abstract>`_.
