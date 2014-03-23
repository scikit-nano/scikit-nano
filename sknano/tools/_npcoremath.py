# -*- coding: utf-8 -*-
"""
==================================================================
Custom NumPy data structures (:mod:`sknano.tools._npcoremath`)
==================================================================

.. currentmodule:: sknano.tools._npcoremath

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np

#from ._corefuncs import check_type

__all__ = ['NPPoint']  # , 'NPVector', 'NPQuaternion']


class NPPoint(np.ndarray):
    """Create a point in :math:`R^3`

    Parameters
    ----------
    x, y, z : float, optional
        :math:`x, y, z` coordinates of point in :math:`R^3` space.
    units : {None, str}, optional
        Units of coordinates.

    """
    def __new__(cls, x=None, y=None, z=None, units=None):
        pass
        #self._p = np.zeros(3, dtype=float)
        #self._units = units
        #for i, pi in enumerate((x, y, z)):
        #    if pi is not None:
        #        self._p[i] = pi
