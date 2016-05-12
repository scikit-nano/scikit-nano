# -*- coding: utf-8 -*-
"""
===============================================================================
Base geometric region classes (:mod:`sknano.core.geometric_regions.base`)
===============================================================================

.. currentmodule:: sknano.core.geometric_regions.base

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod
from functools import reduce, total_ordering

import numpy as np

from sknano.core import BaseClass, TabulateMixin
from sknano.core.math import Point, Points, Vectors, transformation_matrix

# import numpy as np

__all__ = ['GeometricRegion', 'GeometricTransformsMixin']

ndim_errmsg = 'Expected a {}-element array_like object'


class GeometricTransformsMixin:
    """Mixin class providing methods for applying linear algebra \
        transforms to geometric regions."""

    def rotate(self, **kwargs):
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
        # if kwargs.get('axis', None) is None:
        #     kwargs['axis'] = 'z'
        if kwargs.get('anchor_point', None) is None:
            kwargs['anchor_point'] = self.centroid

        if kwargs.get('transform_matrix', None) is None:
            kwargs['transform_matrix'] = transformation_matrix(**kwargs)

        self.points.rotate(**kwargs)
        self.vectors.rotate(**kwargs)

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

    def __deepcopy__(self, memo):
        obj = self.__class__(**self.todict())
        memo[id(self)] = obj
        return obj

    def __eq__(self, other):
        return isinstance(other, type(self)) and \
            self.points == other.points and self.vectors == other.vectors

    def __lt__(self, other):
        return isinstance(other, type(self)) and self.measure < other.measure

    def __str__(self):
        strrep = self._table_title_str()
        objstr = self._obj_mro_str()
        points = self.points
        vectors = self.vectors
        if points:
            title = '.'.join((objstr, points.__class__.__qualname__))
            strrep = '\n'.join((strrep, title, str(points)))
        if vectors:
            title = '.'.join((objstr, vectors.__class__.__qualname__))
            strrep = '\n'.join((strrep, title, str(vectors)))
        return strrep

    @property
    @abstractmethod
    def centroid(self):
        """Centroid of geometric region."""
        raise NotImplementedError

    @property
    @abstractmethod
    def ndim(self):
        """Dimensions of of geometric region."""
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
        """:class:`Point` at minimum extent."""
        # self._update_minmax_points()
        try:
            return self._pmin
        except AttributeError:
            return Point(reduce(lambda x1, x2: np.minimum(x1, x2),
                                self.get_points()))

    @property
    def pmax(self):
        """:class:`Point` at maximum extent."""
        try:
            return self._pmax
        except AttributeError:
            return Point(reduce(lambda x1, x2: np.maximum(x1, x2),
                                self.get_points()))

    def get_points(self):
        """Return list of points from :attr:`GeometricRegion.points` and \
            :attr:`GeometricRegion.vectors`"""
        points = [np.asarray(pt) for pt in self.points]
        points.extend([np.asarray(vec.p) for vec in self.vectors])
        return points

    def center_centroid(self):
        """Center :attr:`~GeometricRegion.centroid` on origin."""
        self.translate(-self.centroid)
