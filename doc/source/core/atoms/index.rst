.. include:: references.txt

.. currentmodule:: sknano.core.atoms

.. _sknano-core-atoms:

===============================================================================
Class representations of nature's building blocks (:mod:`sknano.core.atoms`)
===============================================================================

Introduction
============

The `~sknano.core.atoms` package provides class representations of
atoms and molecules.

Getting Started
================

Of the classes implemented in `~sknano.core.atoms`, there are two base
classes from which most classes in `~sknano.core.atoms` inherit from:
|Atom| and |Atoms|.
The |Atom| class represents a single atom. The |Atoms| class is a container
class for |Atom| class instances. Sub-classes of |Atom| classes
add new atom attributes to the |Atom| class. Every |Atom| sub-class
has a corresponding container class that sub-classes the |Atoms| class.

Base |Atom|/|Atoms| classes
---------------------------
The `Atom` class represents a single atom. The `Atoms` class is a container
class for `Atom` class instances. Sub-classes of `Atom` classes
add new atom attributes to the `Atom` class. Every `Atom` sub-class
has a corresponding container class that sub-classes the `Atoms` class.

Mixins for extending functionality
----------------------------------

The `~sknano.core.atoms.mixins` sub-package provides mixin classes
which extend the existing |Atom| and/or |Atoms| sub-classes.

Classes for molecular dynamics simulations
------------------------------------------

There are two classes for molecular dynamics simulations |Trajectory| and
|Snapshot|.


Using ``atoms``
================

.. toctree::
   :maxdepth: 1

   mixins

Reference/API
=============

.. automodapi:: sknano.core.atoms

.. automodapi:: sknano.core.atoms.atoms

.. automodapi:: sknano.core.atoms.basis_atoms

.. automodapi:: sknano.core.atoms.charged_atoms

.. automodapi:: sknano.core.atoms.cn_atoms

.. automodapi:: sknano.core.atoms.dipole_atoms

.. automodapi:: sknano.core.atoms.energy_atoms

.. automodapi:: sknano.core.atoms.force_atoms

.. automodapi:: sknano.core.atoms.id_atoms

.. automodapi:: sknano.core.atoms.image_atoms

.. automodapi:: sknano.core.atoms.lattice_atoms

.. automodapi:: sknano.core.atoms.md_atoms

.. automodapi:: sknano.core.atoms.molecules

.. automodapi:: sknano.core.atoms.neighbor_atoms

.. automodapi:: sknano.core.atoms.selections

.. automodapi:: sknano.core.atoms.structure_atoms

.. automodapi:: sknano.core.atoms.trajectory

.. automodapi:: sknano.core.atoms.type_atoms

.. automodapi:: sknano.core.atoms.vdW_atoms

.. automodapi:: sknano.core.atoms.velocity_atoms

.. automodapi:: sknano.core.atoms.xyz_atoms
