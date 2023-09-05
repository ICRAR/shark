Changelog
=========

.. rubric:: 2.0.0

* Many changes to the physical models SHARK, which are collectively described in
  Lagos et al. (2023, submitted to MNRAS). There are also many changes in the
  way we create extinction parameters, which are described in the paper Lagos et
  al. (2019, MNRAS, 489.4196).
* Many improvements to the our PSO support.
  This includes better logging of detailed constraint evaluation information,
  improved stability in a few corner cases,
  and offline evaluation of previous shark runs.
* Improved the memory footprint of |s| executions.
  We have made a major overhaul of the code
  to be more memory efficient,
  thus allowing for better resource usage,
  specially in HPC systems
  and PSO executions.
  Several experiments with our mini-SURFS and medi-SURFS datasets
  show a decrease of about 20% on peak memory usage.
* Made |s| fully reproducible in multi-threaded mode,
  even when running a different number of threads
  than a previous execution.

.. rubric:: 1.2.1

* Fixed compilation problem related to OpenMP support.

.. rubric:: 1.2.0

* Initial version of the :doc:`PSO optimization routines <optim>`,
  originally contributed by Mawson Sammons.
  This initial version allows users to run PSO over |s|
  using an arbitrary number of parameters.
  Two evaluation functions (chi2 and log student-t)
  and three constraints with flexible domains (SMF at z=0,1 and HIMF at z=0)
  can be chosen from.
  Execution in both desktop/laptop and HPC environments is supported.
* |s| can also read halo inputs from hydro-dynamical simulations.
* Correctly tracking history of subhalos across snapshots.
* Simple model of radiation pressure QSO feedback included.
* Extended outputs including halo properties.
* Add optional execution parameter to seed random number engines,
  and recording it on output files.
  These two options allow users to fully reproduce a previous |s| run.
  For now this works for single-threaded executions only.
* Improved support for the MSVC compiler.
  |s| now correctly compiles, runs, and standard plots work correctly on Windows.
* Improved the |ss| script to accept additional environment variables
  to set default values for submission parameters.
* Miscelaneous additions and modifications to plotting scripts.

.. rubric:: 1.1.0

* Fixed bug in spin parameter of halos
  as it was randomly chosen despite the value
  the progenitor halo had.
  Now this is coherently propagated through the halo history.
* Type-2 velocity dispersions now consider
  virial theorem at the radius at which the galaxy lives.
* Added new concentration model of Dutton14
* |s| successfully building on Windows with Visual Studio compiler.
* Explicitly supporting OpenMP >= 2.0;
  before we were implicitly requiring >= 3.0.
* Fixed some small problems with |ss|
* Plot scripts simplified to use the redshift table
  to correctly calculate the snapshots they should read

.. rubric:: 1.0.0

* First public release of shark
