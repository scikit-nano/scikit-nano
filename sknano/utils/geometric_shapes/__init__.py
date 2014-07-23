# -*- coding: utf-8 -*-
"""
============================================================================
Abstract geometric data structures (:mod:`sknano.tools.geometric_shapes`)
============================================================================

.. currentmodule:: sknano.tools.geometric_shapes

Contents
========

2D Shapes
----------

.. autosummary::
   :toctree: generated/

   Circle
   Ellipse
   Parallelogram
   Rhombus
   Rhomboid
   Rectangle
   Square
   Polygon

3D Shapes
----------

.. autosummary::
   :toctree: generated/

   Cube
   Cuboid
   Ellipsoid
   Spheroid
   Sphere
   Polyhedron
   Hexahedron
   Parallelepiped
   Rhombohedron

"""
from __future__ import absolute_import, division, print_function

__docformat__ = 'restructuredtext en'

from ._2D_shapes import *
from ._3D_shapes import *

__all__ = [s for s in dir() if not s.startswith('_')]
