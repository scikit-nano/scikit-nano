===============================================================================
Class representations of nature's building blocks (:mod:`sknano.core.atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms

Introduction
============

The `~sknano.core.atoms` package provides class representations of
atoms and molecules.

Getting Started
----------------------
Of the classes implemented in `~sknano.core.atoms`, there are two base
classes from which most classes in `~sknano.core.atoms` inherit from:
`~sknano.core.atoms.Atom` and `~sknano.core.atoms.Atoms`.
The `Atom` class represents a single atom. The `Atoms` class is a container
class for `Atom` class instances. Sub-classes of `Atom` classes
add new atom attributes to the `Atom` class. Every `Atom` sub-class
has a corresponding container class that sub-classes the `Atoms` class.


Mixins for extending functionality
----------------------------------
    * `~sknano.core.atoms.mixins`

Classes for molecular dynamics simulations
------------------------------------------

There are two classes for molecular dynamics simulations
`~sknano.core.atoms.Trajectory` and `~sknano.core.atoms.Snapshot`
classes

Helper functions for atom objects
---------------------------------

   * `~sknano.core.atoms.compute_angle`
   * `~sknano.core.atoms.compute_bond`
   * `~sknano.core.atoms.compute_dihedral`
   * `~sknano.core.atoms.compute_improper`
   * `~sknano.core.atoms.vdw_radius_from_basis`

Reference/API
=============

.. automodapi:: sknano.core.atoms
