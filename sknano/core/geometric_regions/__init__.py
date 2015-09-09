# -*- coding: utf-8 -*-
"""
============================================================================
Abstract geometric data structures (:mod:`sknano.core.geometric_regions`)
============================================================================

.. currentmodule:: sknano.core.geometric_regions

Contents
========

Base/mixin classes for all geometric regions
-----------------------------------------------------------
.. autosummary::
   :toctree: generated/

   GeometricRegion
   Geometric2DRegion
   Geometric3DRegion
   GeometricTransformsMixin

2D Regions
----------
.. autosummary::
   :toctree: generated/

   Parallelogram
   Rectangle
   Square
   Ellipse
   Circle
   Triangle

3D Regions
----------

.. autosummary::
   :toctree: generated/

   Parallelepiped
   Cuboid
   Cube
   Ellipsoid
   Sphere
   Cylinder
   Cone

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from ._base import *
from ._funcs import *
from ._2D_regions import *
from ._3D_regions import *

__all__ = [s for s in dir() if not s.startswith('_')]
