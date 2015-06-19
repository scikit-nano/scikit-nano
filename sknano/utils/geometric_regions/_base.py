# -*- coding: utf-8 -*-
"""
===============================================================================
Base classes for geometric regions (:mod:`sknano.utils.geometric_shapes._base`)
===============================================================================

.. currentmodule:: sknano.utils.geometric_shapes._base

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod
from sknano.core.math import Points, Vectors, transformation_matrix

import numpy as np

__all__ = ['GeometricRegion', 'GeometricTransformsMixin']


class GeometricTransformsMixin:

    def rotate(self, angle=None, axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None,
               degrees=False, transform_matrix=None, verbose=False, **kwargs):
        """Rotate `GeometricRegion` `Points` and `Vectors`.

        Parameters
        ----------
        angle : float
        axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        degrees : bool, optional
        transform_matrix : :class:`~numpy:numpy.ndarray`

        """
        if transform_matrix is None:
            transform_matrix = \
                transformation_matrix(angle=angle, axis=axis,
                                      anchor_point=anchor_point,
                                      rot_point=rot_point,
                                      from_vector=from_vector,
                                      to_vector=to_vector, degrees=degrees,
                                      verbose=verbose, **kwargs)

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


class GeometricRegion(metaclass=ABCMeta):
    """Abstract base class for geometric regions."""

    def __init__(self):
        self.points = Points()
        self.vectors = Vectors()
        self.fmtstr = ""

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    @property
    @abstractmethod
    def centroid(self):
        """Centroid of geometric region."""
        raise NotImplementedError

    @abstractmethod
    def contains(self):
        """Check if point is contained within geometric region."""
        raise NotImplementedError

    @abstractmethod
    def todict(self):
        """Return `dict` of `__init__` parameters."""
        raise NotImplementedError
