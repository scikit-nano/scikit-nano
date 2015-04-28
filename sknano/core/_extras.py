# -*- coding: utf-8 -*-
"""
===============================================================================
Misc core functions, constants, etc. (:mod:`sknano.core._extras`)
===============================================================================

.. currentmodule:: sknano.core._extras

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

__docformat__ = 'restructuredtext en'

import numpy as np

__all__ = ['components', 'dimensions', 'xyz', 'xyz_axes',
           'dedupe', 'rezero_array']

components = dimensions = xyz = xyz_axes = ('x', 'y', 'z')


def dedupe(items, key=None):
    """Remove duplicate values in a sequence, but preserve order of remaining \
        items.

    Parameters
    ----------
    items : sequence
    key : {None, function}, optional
        function that converts sequence items into a hashable type for
        the purposes of duplicate detection.

    Returns
    -------
    items : set

    Examples
    --------
    >>> a = [{'x': 1, 'y': 2}, {'x': 1, 'y': 3},
    ...      {'x': 1, 'y': 2}, {'x': 2, 'y': 4}]
    >>> list(dedupe(a, key=lambda d: (d['x'], d['y'])))
    [{'x': 1, 'y': 2}, {'x': 1, 'y': 3}, {'x': 2, 'y': 4}]
    >>> list(dedupe(a, key=lambda d: d['x']))
    [{'x': 1, 'y': 2}, {'x': 2, 'y': 4}]

    """
    seen = set()
    for item in items:
        val = item if key is None else key(item)
        if val not in seen:
            yield item
            seen.add(val)


def rezero_array(a, epsilon=5.0*np.finfo(float).eps):
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
