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
   :toctree: api/

   Point
   Points
   Vector
   Vectors

Linear algebra transforms
----------------------------
.. autosummary::
   :toctree: api/

   rotate
   Rx
   Ry
   Rz
   reflection_matrix
   rotation_matrix
   scaling_matrix
   transformation_matrix
   axis_angle_from_rotation_matrix
   rotation_matrix_from_augmented_matrix
   translation_from_augmented_matrix

Number theory
---------------
.. autosummary::
   :toctree: api/

   totient_func

"""
from __future__ import absolute_import
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from .extras import *
from .point import *
from .points import *
from .vector import *
from .vectors import *
from .quaternion import *
from .transforms import *

__all__ = [s for s in dir() if not s.startswith('_')]

from . import point as point
from . import vector as vector
from . import transforms as transforms

__all__ += ['point', 'vector', 'transforms']
