Output files
============

For each snapshot of an execution,
|s| optionally outputs three different types of files:

 * ``galaxies.hdf5`` contains the information of galaxies
 * ``star_formation_histories.hdf5`` contains the star formation history of the galaxies

We list here each of the HDF5 groups and datasets
found on each of these files.
This list is automatically calculated from the files themselves.

.. _output.galaxies:

Galaxies
--------

.. include:: hdf5_properties/galaxies.rst


Star formation histories
------------------------

.. include:: hdf5_properties/star_formation_histories.rst
