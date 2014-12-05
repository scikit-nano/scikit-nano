# -*- coding: utf-8 -*-
"""
===============================================================================
Misc core functions, constants, etc. (:mod:`sknano.core._extras`)
===============================================================================

.. currentmodule:: sknano.core._extras

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

from itertools import chain, tee

from six.moves import zip

__docformat__ = 'restructuredtext en'

import numpy as np

__all__ = ['components', 'dimensions', 'xyz', 'xyz_axes',
           'AttrDict', 'cyclic_pairs', 'rezero_array']

components = dimensions = xyz = xyz_axes = ('x', 'y', 'z')


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def cyclic_pairs(iterable):
    """Generate a cyclic list of all pair subsequences of elements from the \
        input `iterable`.

    Parameters
    ----------
    iterable : sequence

    Returns
    -------
    :class:`~python:list`

    """
    a, b = tee(iterable)
    return list(zip(a, chain(b, [next(b)])))


def rezero_array(a, epsilon=None):
    """Rezero elements of array `a` with absolute value \
        *less than or equal to* `epsilon`.

    Parameters
    ----------
    a : :class:`~numpy:numpy.ndarray`
    epsilon : {None, float}, optional

    Returns
    -------
    a : :class:`~numpy:numpy.ndarray`

    """
    if not isinstance(a, np.ndarray):
        raise TypeError('Expected a numpy array')

    if epsilon is None:
        epsilon = np.finfo(float).eps

    a[np.where(np.abs(a) <= epsilon)] = 0.0

    return a
