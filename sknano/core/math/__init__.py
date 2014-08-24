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
   Vector

Linear algebra transforms
----------------------------
.. autosummary::
   :toctree: generated/

   rotation_matrix
   transformation_matrix
   rotation_transform

Number theory
---------------
.. autosummary::
   :toctree: generated/

   totient_func

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._funcs import *
from ._transforms import *

from ._point import *
from ._vector import *

__all__ = [s for s in dir() if not s.startswith('_')]

from . import _point as point
from . import _vector as vector

__all__ += ['point', 'vector']
