# -*- coding: utf-8 -*-
"""
===============================================================================
Base classes for geometric regions (:mod:`sknano.utils.geometric_shapes._base`)
===============================================================================

.. currentmodule:: sknano.utils.geometric_shapes._base

"""
from __future__ import absolute_import, division, print_function
import six
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod, abstractproperty
from sknano.core.math import Points, Vectors, transformation_matrix

import numpy as np

__all__ = ['GeometricRegion', 'GeometricTransformsMixin']


class GeometricTransformsMixin(object):

    def rotate(self, angle=None, rot_axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None,
               deg2rad=False, transform_matrix=None, verbose=False):
        """Rotate `GeometricRegion` `Points` and `Vectors`.

        Parameters
        ----------
        angle : float
        rot_axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        deg2rad : bool, optional
        transform_matrix : :class:`~numpy:numpy.ndarray`

        """
        if transform_matrix is None:
            transform_matrix = \
                transformation_matrix(angle=angle, rot_axis=rot_axis,
                                      anchor_point=anchor_point,
                                      rot_point=rot_point,
                                      from_vector=from_vector,
                                      to_vector=to_vector, deg2rad=deg2rad,
                                      verbose=verbose)

        self.points.rotate(transform_matrix=transform_matrix)
        self.vectors.rotate(transform_matrix=transform_matrix)

    def translate(self, t, fix_anchor_points=False):
        """Translate `GeometricRegion`\ s `Points` and `Vectors` by \
            :class:`~sknano.core.math.Vector` `t`.

        Parameters
        ----------
        t : :class:`~sknano.core.math.Vector`
        fix_anchor_points : bool, optional

        """
        self.points.translate(t)
        self.vectors.translate(t, fix_anchor_points=fix_anchor_points)


class GeometricRegion(six.with_metaclass(ABCMeta, object)):
    """Abstract base class for geometric regions."""

    def __init__(self):
        self.points = Points()
        self.vectors = Vectors()
        self.limits = {'x': {'min': -np.inf, 'max': np.inf},
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
