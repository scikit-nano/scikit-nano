# -*- coding: utf-8 -*-
"""
============================================================================
Tools for testing (:mod:`sknano.testing`)
============================================================================

.. currentmodule:: sknano.testing

Contents
========

Helper funcs to generate test data
---------------------------------------

.. autosummary::
   :toctree: generated/

   generate_atoms
   generate_structure

Tools for timing processes
---------------------------

.. autosummary::
   :toctree: generated/

   Timer

Test fixtures
--------------

.. autosummary::
   :toctree: generated/

   AtomsTestFixture
   GeneratorTestFixture
   IOTestFixture
   TempfileTestFixture
   DUMPTestFixture
   GeometricRegionsTestFixture
   Geometric2DRegionsTestFixture
   Geometric3DRegionsTestFixture

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from ._funcs import generate_atoms, generate_structure
from ._timer import Timer
from ._tools import AtomsTestFixture, GeneratorTestFixture, IOTestFixture, \
    TempfileTestFixture, DUMPTestFixture, GeometricRegionsTestFixture, \
    Geometric2DRegionsTestFixture, Geometric3DRegionsTestFixture

__all__ = [s for s in dir() if not s.startswith('_')]
