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
* :class:`~sknano.generators.SWNTBundleGenerator`
* :class:`~sknano.generators.MWNTGenerator`
* :class:`~sknano.generators.MWNTBundleGenerator`


.. ipython::

   In [1]: from sknano.generators import SWNTBundleGenerator

   In [2]: bundle = SWNTBundleGenerator(n=10, m=5, nz=2, bundle_geometry='hexagon')

   In [3]: bundle.save_data()
