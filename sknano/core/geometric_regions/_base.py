# -*- coding: utf-8 -*-
"""
===============================================================================
Base geometric region classes (:mod:`sknano.core.geometric_regions._base`)
===============================================================================

.. currentmodule:: sknano.core.geometric_regions._base

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod
from functools import reduce, total_ordering

import numpy as np

from sknano.core import BaseClass, TabulateMixin
from sknano.core.math import Points, Vectors, transformation_matrix

# import numpy as np

__all__ = ['GeometricRegion', 'GeometricTransformsMixin']


class GeometricTransformsMixin:
    """Mixin class providing methods for applying linear algebra \
        transforms to geometric regions."""

    def rotate(self, angle=None, axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None,
               degrees=False, transform_matrix=None, verbose=False, **kwargs):
        """Rotate `GeometricRegion` :attr:`~GeometricRegion.points` and \
            :attr:`~GeometricRegion.vectors`.

        Parameters
        ----------
        angle : float
        axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        degrees : bool, optional
        transform_matrix : :class:`~numpy:numpy.ndarray`


        See Also
        --------
        sknano.core.math.rotate

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
        """Translate `GeometricRegion` :attr:`~GeometricRegion.points` and \
            :attr:`~GeometricRegion.vectors` by \
            :class:`~sknano.core.math.Vector` `t`.

        Parameters
        ----------
        t : :class:`~sknano.core.math.Vector`
        fix_anchor_points : bool, optional

        See Also
        --------
        sknano.core.math.translate

        """
        self.points.translate(t)
        self.vectors.translate(t, fix_anchor_points=fix_anchor_points)


@total_ordering
class GeometricRegion(BaseClass, TabulateMixin, metaclass=ABCMeta):
    """Abstract base class for all geometric regions.

    Attributes
    ----------
    points : :class:`Points`
        Collection of all :class:`~sknano.core.math.Point` objects
        defining the :class:`GeometricRegion`.
    vectors : :class:`Vectors`
        Collection of all :class:`~sknano.core.math.Vector` objects
        defining the :class:`GeometricRegion`

    """
    def __init__(self):
        super().__init__()
        self.points = Points()
        self.vectors = Vectors()

    def __eq__(self, other):
        return isinstance(other, type(self)) and \
            self.points == other.points and self.vectors == other.vectors

    def __lt__(self, other):
        return isinstance(other, type(self)) and self.measure < other.measure

    def __str__(self):
        strrep = self._table_title_str()
        points = self.points
        vectors = self.vectors
        if points:
            strrep = '\n'.join((strrep, str(points)))
        if vectors:
            strrep = '\n'.join((strrep, str(vectors)))
        return strrep

    @property
    @abstractmethod
    def centroid(self):
        """Centroid of geometric region."""
        raise NotImplementedError

    @property
    def center(self):
        """Alias for :attr:`~GeometricRegion.centroid`."""
        return self.centroid

    @property
    @abstractmethod
    def measure(self):
        """Measure of geometric region."""
        raise NotImplementedError

    @abstractmethod
    def contains(self, point):
        """Test region membership of `point` in :class:`GeometricRegion`."""
        raise NotImplementedError

    @property
    def pmin(self):
        """:class:`Point` at minimum extent in :math:`(x,y,z)` dimensions."""
        points = [pt for pt in self.points]
        points.extend([vec.p for vec in self.vectors])
        return reduce(lambda x1, x2: np.minimum(x1, x2), points)

    @property
    def pmax(self):
        """:class:`Point` at maximum extent in :math:`(x,y,z)` dimensions."""
        points = [pt for pt in self.points]
        points.extend([vec.p for vec in self.vectors])
        return reduce(lambda x1, x2: np.maximum(x1, x2), points)

    def center_centroid(self):
        """Center :attr:`~GeometricRegion.centroid` on origin."""
        self.translate(-self.centroid)
