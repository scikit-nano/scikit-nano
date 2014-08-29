# -*- coding: utf-8 -*-
"""
===============================================================================
Base classes for geometric regions (:mod:`sknano.utils.geometric_shapes._base`)
===============================================================================

.. currentmodule:: sknano.utils.geometric_shapes._base

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod, abstractproperty
#from sknano.core.math import rotation_matrix

import numpy as np

__all__ = ['GeometricRegion', 'GeometricTransformsMixin']


class GeometricTransformsMixin(object):

    def rotate(self, angle, rot_axis=None, anchor_point=None, deg2rad=False):
        """Rotate region about centroid or arbitrary vector."""
        for p in self.points:
            p.rotate(angle, rot_axis=rot_axis, anchor_point=anchor_point,
                     deg2rad=deg2rad)
        for v in self.vectors:
            v.rotate(angle, rot_axis=rot_axis, anchor_point=anchor_point,
                     deg2rad=deg2rad)

    def translate(self, t):
        """Translate region."""
        for p in self.points:
            p.translate(t)

        for v in self.vectors:
            v.translate(t)


class GeometricRegion(object):
    """Abstract base class for geometric regions."""
    __metaclass__ = ABCMeta

    def __init__(self):
        self._points = []
        self._vectors = []
        self._limits = {'x': {'min': -np.inf, 'max': np.inf},
                        'y': {'min': -np.inf, 'max': np.inf},
                        'z': {'min': -np.inf, 'max': np.inf}}

    @abstractproperty
    def centroid(self):
        """Centroid of geometric region."""
        raise NotImplementedError

    @abstractmethod
    def contains_point(self):
        """Check if point is contained within geometric region."""
        raise NotImplementedError

    @property
    def limits(self):
        return self._limits

    @property
    def points(self):
        return self._points

    @property
    def vectors(self):
        return self._vectors
