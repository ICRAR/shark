Configuration
=============

.. _config.file:

Via a configuration file
------------------------

|s| uses a configuration file as its primary source of configuration.
Configuration files are divided in sections,
each containing a number of configuration items.
Lines starting with a ``#`` are considered comments,
and can appear anywhere in the file.
Blank lines are ignored.

Sections start with a line like ``[abc]``,
where ``abc`` is the name of the section.
After that zero or more configuration items can appear,
each in the form of ``name = value``.
Spaces before or after the ``=`` sign are optional,
but encouraged for readability.

For example, this is a small sample configuration file::

 # Sample configuration file
 #
 # General execution parameters
 [execution]
 output_format = hdf5
 output_snapshots = 199-55

 # Cosmology parameters
 [cosmology]
 omega_m = 0.3121

|s| allows more than one configuration file
to be used at the same time.
If options appear in more than one configuration file,
the last configuration file specified in the command line
takes precedence.

.. _config.cmdline:

Via the command line
--------------------

Additionally, |s| allows users to pass configuration values
via command line options.
Configuration values are given through the command line
using the ``-o`` option, like this::

 $> shark -o cosmology.omega_m=0.4 ...

Note how configuration item name has the form ``<group>.<name>``.
More than one ``-o`` option can be given to pass down
more than one configuration value.
For further information on how to run |s|
see :doc:`running`.

Configuration values given on the command line
take precedence over values specified
on any configuration file.
In the example above,
and if using the sample configuration shown in `config.file`,
|s| will effectively run with a value of ``0.4``
for ``cosmology.omega_m``.
