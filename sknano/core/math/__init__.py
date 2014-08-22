# -*- coding: utf-8 -*-
"""
======================================================================
Core package modules (:mod:`sknano.core`)
======================================================================

.. currentmodule:: sknano.core

Contents
========

Abstract data structures for math
--------------------------------------
.. autosummary::
   :toctree: generated/

   Point
   Vector

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._point import *
from ._vector import *

__all__ = [s for s in dir() if not s.startswith('_')]

from . import _point as point
from . import _vector as vector

__all__ += ['point', 'vector']
