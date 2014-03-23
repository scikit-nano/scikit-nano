# -*- coding: utf-8 -*-
"""
===============================================================================
Linear algebra functions for transformations (:mod:`sknano.tools._transforms`)
===============================================================================

.. currentmodule:: sknano.tools._transforms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np

__all__ = ['rotation_matrix']


def rotation_matrix(angle, rot_axis='z', deg2rad=False):
    """Generate a rotation matrix.

    Parameters
    ----------
    angle : float
        rotation angle in radians
    rot_axis : {'x', 'y', 'z'}, optional
        rotation axis
    deg2rad : bool, optional
        Angle is in degrees and needs to be converted to radians

    Returns
    -------
    ndarray
        rotation matrix

    """
    if deg2rad:
        angle = np.radians(angle)
    if rot_axis == 'x':
        return np.array([[1, 0, 0],
                         [0, np.cos(angle), -np.sin(angle)],
                         [0, np.sin(angle), np.cos(angle)]])
    elif rot_axis == 'y':
        return np.array([[np.cos(angle), 0, np.sin(angle)],
                         [0, 1, 0],
                         [-np.sin(angle), 0, np.cos(angle)]])
    else:
        return np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0],
                         [0, 0, 1]])
