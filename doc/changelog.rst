Changelog
=========

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
