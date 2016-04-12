.. include:: references.txt

.. currentmodule:: sknano.generators

.. _sknano-generators:

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
|CrystalStructureGenerator|, while
nanostructure generator classes inherit from
|NanoStructureGenerator| base class.

Both |CrystalStructureGenerator| and
|NanoStructureGenerator| inherit from
|GeneratorBase| base class.

The |GeneratorMixin| class provides a simple concrete
implementation of the `~sknano.generators.GeneratorBase.generate` method.


Example::

    >>> from sknano.generators import SWNTGenerator

Reference/API
=============

.. automodapi:: sknano.generators
