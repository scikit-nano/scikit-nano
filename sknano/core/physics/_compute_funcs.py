# -*- coding: utf-8 -*-
"""
===============================================================================
Physics functions (:mod:`sknano.core.physics._compute_funcs`)
===============================================================================

.. currentmodule:: sknano.core.physics._compute_funcs

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core.math import Vector

__all__ = ['compute_centroid', 'compute_inertia_tensor']


def compute_centroid(coords):
    """Compute the centroid of coordinates.

    .. math::
       \\mathbf{C} =
       \\frac{\\sum_{i=1}^{N}m_i\\mathbf{r}_i}{\\sum_{i=1}^{N}m_i}

    Returns
    -------
    C : `~sknano.core.math.Vector`
        The position vector of the centroid coordinates.
    """
    C = Vector(np.mean(coords, axis=0))
    C.rezero()
    return C


def compute_inertia_tensor(masses, coords):
    """Compute the inertia tensor."""
    x, y, z = [coords[:, i] for i in range(3)]

    Ixx = (masses * (y**2 + z**2)).sum()
    Iyy = (masses * (x**2 + z**2)).sum()
    Izz = (masses * (x**2 + y**2)).sum()
    Ixy = Iyx = (-masses * x * y).sum()
    Ixz = Izx = (-masses * x * z).sum()
    Iyz = Izy = (-masses * y * z).sum()
    return np.array([[Ixx, Ixy, Ixz], [Iyx, Iyy, Iyz], [Izx, Izy, Izz]])
