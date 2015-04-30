# -*- coding: utf-8 -*-
"""
====================================================================
Custom NumPy Quaternion class (:mod:`sknano.core.math._quaternion`)
====================================================================

.. currentmodule:: sknano.core.math._quaternion

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import copy
import numbers
import warnings

import numpy as np
np.seterr(all='warn')

# from ._point import Point
# from ._transforms import rotate, transformation_matrix

__all__ = ['Quaternion']


class Quaternion(np.ndarray):
    """Abstract object representation of a quaternion.

    Parameters
    ----------
    dtype : data-type, optional
    copy : bool, optional

    Examples
    --------

    """
    __array_priority__ = 15.0
    _verbosity = 0

    def __new__(cls, dtype=None, copy=True):
        pass
