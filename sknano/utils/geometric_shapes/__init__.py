# -*- coding: utf-8 -*-
"""
============================================================================
Abstract geometric data structures (:mod:`sknano.utils.geometric_shapes`)
============================================================================

.. currentmodule:: sknano.utils.geometric_shapes

Contents
========

Abstract Base Classes for defining new geometric shapes
--------------------------------------------------------
.. autosummary::
   :toctree: generated/

   GeometricRegion
   Geometric2DRegion
   Geometric3DRegion

2D Shapes
----------
.. autosummary::
   :toctree: generated/

   Parallelogram
   Rectangle
   Square
   Ellipse
   Circle

3D Shapes
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

__docformat__ = 'restructuredtext en'

from ._base import *
from ._2D_shapes import *
from ._3D_shapes import *

__all__ = [s for s in dir() if not s.startswith('_')]
