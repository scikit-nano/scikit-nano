# -*- coding: utf-8 -*-
"""
============================================================================
Testing tools (:mod:`sknano.testing`)
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
   GeneratorTestFixtures

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from ._funcs import *
from ._timer import *
from ._tools import *

__all__ = [s for s in dir() if not s.startswith('_')]
