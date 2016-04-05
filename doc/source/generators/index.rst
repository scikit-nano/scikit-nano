.. include:: references.txt

.. _sknano-generators:

.. currentmodule:: sknano.generators

======================================================================
Structure generators (:mod:`sknano.generators`)
======================================================================

Introduction
=============

The `sknano.generators` package provides classes for generating structure data
for crystal structures and nanostructures.

Getting Started
===============

Crystal structure generator classes are sub-classes of
`~sknano.generators.CrystalStructureGenerator`, while
nanostructure generator classes inherit from
`~sknano.generators.NanoStructureGenerator` base class.

Both `~sknano.generators.CrystalStructureGenerator` and
`~sknano.generators.NanoStructureGenerator` inherit from
`~sknano.generators.GeneratorBase` base class.

The `~sknano.generators.GeneratorMixin` class provides a simple concrete
implementation of the `~sknano.generators.GeneratorBase.generate` method.


Reference/API
=============

.. automodapi:: sknano.generators
