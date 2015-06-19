# -*- coding: utf-8 -*-
"""
============================================================================
Abstract geometric data structures (:mod:`sknano.core.geometric_regions`)
============================================================================

.. currentmodule:: sknano.core.geometric_regions

Contents
========

Abstract Base Classes for defining new geometric regions
--------------------------------------------------------
.. autosummary::
   :toctree: generated/

   GeometricRegion
   Geometric2DRegion
   Geometric3DRegion

2D Regions
----------
.. autosummary::
   :toctree: generated/

   Parallelogram
   Rectangle
   Square
   Ellipse
   Circle

3D Regions
----------

.. autosummary::
   :toctree: generated/

   Parallelepiped
   Cuboid
   Cube
   Ellipsoid
   Spheroid
   Sphere

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from ._base import *
from ._2D_regions import *
from ._3D_regions import *

__all__ = [s for s in dir() if not s.startswith('_')]
