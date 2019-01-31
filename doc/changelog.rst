Changelog
=========

.. rubric:: Development version

* Add optional execution parameter to seed random number engines,
  and recording it on output files.
  These two options allow users to fully reproduce a previous |s| run.
  For now this works for single-threaded executions only.
* Improved support for the MSVC compiler.
  |s| now correctly compiles, runs, and standard plots work correctly on Windows.
* Improved the |ss| script to accept additional environment variables
  to set default values for submission parameters.

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
