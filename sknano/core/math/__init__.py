# -*- coding: utf-8 -*-
"""
======================================================================
Core data structures for math (:mod:`sknano.core.math`)
======================================================================

.. currentmodule:: sknano.core.math

Contents
========

Abstract object representations for points and vectors
-------------------------------------------------------
.. autosummary::
   :toctree: generated/

   Point
   Points
   Vector
   Vectors

Linear algebra transforms
----------------------------
.. autosummary::
   :toctree: generated/

   rotate
   Rx
   Ry
   Rz
   reflection_matrix
   rotation_matrix
   scaling_matrix
   transformation_matrix
   axis_angle_from_rotation_matrix

Number theory
---------------
.. autosummary::
   :toctree: generated/

   totient_func

"""
from __future__ import absolute_import
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._extras import *
from ._point import *
from ._points import *
from ._vector import *
from ._vectors import *
from ._quaternion import *
from ._transforms import *

__all__ = [s for s in dir() if not s.startswith('_')]

from . import _point as point
from . import _vector as vector
from . import _transforms as transforms

__all__ += ['point', 'vector', 'transforms']
