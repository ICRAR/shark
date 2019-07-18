Constraints implementation
==========================

Constraints in |s|'s PSO implementation
are currently independent from each other,
so they can be selected individually
and even weighted separately.
This makes adding new constraints *relatively* easy,
as only the constraint itself has to be implemented,
while the rest of the code doesn't require further changes.

Constraints are implemented in the ``constraints.py`` module.
Each one is implemented as a new class
deriving from the ``Constraint`` class,
which provides the common functionality for all of them
(e.g., plotting, domain filtering, interpolation, etc).
Subclasses only need to take care of loading
the relevant observational and model data,
and provide it to the parent class.
Subclasses also need to inform the parent class
what is their default domain,
although this can be changed explicitly later
by the user.

In summary, three things are required from subclasses:

 * A default domain specification
 * The observational data
 * The model data


.. _constraints.new:

Implementing a new constraint
-----------------------------

To implement a new constraint,
a new class must be written in the ``constraints.py`` module
that derives from the ``Constraint`` class
and that offers the data required by the parent class.
Here is an example of a minimal class:

.. code-block:: python

 class MyConstraint(Constraint):

     domain = [3, 10]
     z = [0]

     def get_obs_x_y_err(self, h0):
         x = [4, 5.6, 7]
         y = [10, 11, 12]
         err_dn = [0.1, 0.2, 0.1]
         err_up = [0.2, 0.34, 0.2]
         return x, y, err_dn, err_up

     def get_model_x_y(self, hist_smf, hist_HImf)
         x = [4, 6, 8]
         y = [12, 11, 13]
         return x, y

Please note how in this example
the observational and model data
do not fall necessarily under the same domain
(e.g, there is a model data point at ``x = 6``,
but no corresponding point in the observational data for it).
This difference in domain
should not be a worry at this level;
it is dealt with later on the parent ``Constraint`` class,
and therefore the implementation of these methods
should be as simple as possible.


Specifying a domain
^^^^^^^^^^^^^^^^^^^

The easiest bit is the default domain specification.
This is specified as a two-element list
with the lower bound first, and the upper bound second.
The actual meaning of these numbers (i.e., their units)
will depend on the data being handled by the constraint.
This list must be called ``domain``
and must be specified at the class level.

Loading observations
^^^^^^^^^^^^^^^^^^^^

Loading observational data is done in the ``get_obs_x_y_err`` method.
This method must return four elements, each an array,
with the values for the domain, observations
and the upper and lower errors associated to each observation.

To easily load an observation
from one of the datasets shipped with |s|
(i.e., inside the ``data`` directory)
one can use the ``load_observation`` method
of the parent ``Constraint`` class.
For examples see the current constraint classes,
all of which load data using this mechanism.

The whole data set relevant for the constraint should be returned;
in particular no domain filtering should (or even can)
be performed at this stage.
Error values must be **relative** to the observational values,
not absolute.

Some observational data needs to be adjusted
for the cosmology of the model, for which the value of ``h0``
is provided as an argument.

Loading model data
^^^^^^^^^^^^^^^^^^

Loading model data is done in the ``get_model_x_y`` method.
This is the hardest bit of all,
as it's the more convoluted in the code at the moment.
Loading model data involves:

 * Reading data from one or more HDF5 files
   (i.e., shark's output)
   for one or more redshifts.
 * Manipulating this data to get the final required result.

While ideally this should happen all in isolation
within each individual ``get_model_x_y`` method,
this process is currently broken down
into several parts of the code.
What is currently done
is that the main ``Constraint`` class
is effectively reading and handing over
**all** the data required by **all** constraints
to **all** constraints via the ``get_model_x_y`` method.
Individual ``Constraint`` subclasses then only *choose*
which bits of the overall dataset
they need to read and return.

Moreover, and in order to be able to reuse
some of the logic in the ``standard_plots`` package
(in particular the ``smf`` module),
the code in the ``Constraint`` class loading the model data
looks weird to the eye.

Model data is also most likely dependent
on a particular redshift,
while |s| outputs are organized by snapshot.
The conversion from redshift to the corresponding snapshot
for the given model
can be done via the ``self.redshift_table`` object
(i.e. ``snapshot = self.redshift_table[z]``.

All of the above is the current state of things,
and doesn't mean that it cannot be changed
-- only that some effort is involved in that.
If you are thus looking to implement a new constraint,
please be aware of this added complexity.


Registering new constraint
--------------------------

Finally, the constraint has to be registered
into the list of all known constraints.
This is currently done by adding an entry
into the ``_constraints`` dictionary
of the ``parse`` method.
Simply add a new entry with the name you wish to give your constraint,
followed by the class name implementing it.
Following the example above,
it would look something like this:

.. code-block:: python

 _constraints = {
     'HIMF': HIMF,
     'SMF_z0': SMF_z0,
     'SMF_z1': SMF_z1,
     'MY_CONSTRAINT': MyConstraint,
 }

 After that,
 it should be possible to use the new constraint
 either when running an actual PSO execution,
 or when performing offline evaluation of model outputs.
