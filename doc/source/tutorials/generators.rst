.. _generators_tutorial:

==============================
Structure generators tutorial
==============================

.. sectionauthor:: Andrew Merrill <androomerrill@gmail.com>

The :mod:`~sknano.generators` module provides the following
classes for generating nanostructures:

Fullerene Structure Generators
=================================
* :class:`~sknano.generators.FullereneGenerator`

Graphene Structure Generators
=================================
* :class:`~sknano.generators.GrapheneGenerator`
* :class:`~sknano.generators.BilayerGrapheneGenerator`
* :class:`~sknano.generators.UnrolledSWNTGenerator`

Nanotube Structure Generators
=================================
* :class:`~sknano.generators.SWNTGenerator`
* :class:`~sknano.generators.MWNTGenerator`

For example, if you have a list of chiralities that you want structure data
for, you can do something like this from within an interactive session::

    >>> from sknano.core.structures import generate_Ch_list
    >>> from sknano.generators import SWNTGenerator
    >>> # Generate your list of (n, m) chiralities
    >>> Ch_list = generate_Ch_list(ni=5, nf=25, mi=0, mf=25, handedness='right')
    >>> for Ch in Ch_list:
    ...     SWNTGenerator(Ch).save(structure_format='data')


.. ipython::

   In [1]: from sknano.generators import SWNTGenerator

   In [2]: bundle = SWNTGenerator((10, 5), nz=2, bundle_geometry='hexagon')

   In [3]: bundle.save()

Bulk Structure Generators
===========================
* :class:`~sknano.generators.AlphaQuartzGenerator`
* :class:`~sknano.generators.DiamondGenerator`
* :class:`~sknano.generators.CsClGenerator`
* :class:`~sknano.generators.NaClGenerator`
* :class:`~sknano.generators.ZincblendeGenerator`
* :class:`~sknano.generators.BCCGenerator`
* :class:`~sknano.generators.FCCGenerator`
* :class:`~sknano.generators.MoS2Generator`

Composite Structure Generators
===============================
* :class:`~sknano.generators.LayeredStructureGenerator`
