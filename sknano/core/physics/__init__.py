# -*- coding: utf-8 -*-
"""
======================================================================
Core modules for physics (:mod:`sknano.core.physics`)
======================================================================

.. currentmodule:: sknano.core.physics

Contents
========

Rigid body calculations
-----------------------
.. autosummary::
   :toctree: generated/

   compute_centroid
   compute_inertia_tensor

"""
from __future__ import absolute_import
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._compute_funcs import *

__all__ = [s for s in dir() if not s.startswith('_')]
