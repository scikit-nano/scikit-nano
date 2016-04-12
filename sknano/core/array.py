# -*- coding: utf-8 -*-
"""
===============================================================================
Helper funcs for arrays (:mod:`sknano.core.array`)
===============================================================================

.. currentmodule:: sknano.core.array

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

__docformat__ = 'restructuredtext en'

import numpy as np

__all__ = ['minmax', 'rezero', 'rezero_array']


def minmax(a):
    """Return :class:`~python:tuple` of (min, max) values in `a`."""
    return np.min(a), np.max(a)


def rezero(a, epsilon=5.0*np.finfo(float).eps):
    """Rezero elements of array `a` with absolute value \
        *less than or equal to* `epsilon`.

    Parameters
    ----------
    a : :class:`~numpy:numpy.ndarray`
    epsilon : float, optional

    Returns
    -------
    a : :class:`~numpy:numpy.ndarray`

    """
    if not isinstance(a, np.ndarray):
        raise TypeError('Expected a numpy array')

    a[np.where(np.abs(a) <= epsilon)] = 0.0

    return a

rezero_array = rezero
