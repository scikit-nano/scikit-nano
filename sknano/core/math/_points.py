# -*- coding: utf-8 -*-
"""
==============================================================================
Points class (:mod:`sknano.core.math._points`)
==============================================================================

.. currentmodule:: sknano.core.math._points

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from sknano.core import UserList, TabulateMixin
from ._transforms import transformation_matrix
# from sknano.core.geometric_regions import Cuboid  # , Rectangle
from ._point import Point

__all__ = ['Points']

operand_shape_error_msg = \
    "operands could not be broadcast together with shapes {}, {}"


class Points(TabulateMixin, UserList):
    """Container class for collection of `Point` objects.

    Parameters
    ----------
    points : {None, sequence, `Points`}, optional
        if not `None`, then a list of `Point` instance objects or an
        existing `Points` instance object.

    """
    def __init__(self, points=None):
        super().__init__(initlist=points)
        self.fmtstr = "{points!r}"

    def _tabular_data(self):
        begin = len('Point(')
        fmt = super()._tabular_data_format_string
        values = list(zip(['P{}'.format(i+1) for i in range(len(self))],
                          [fmt(pt, begin, end=-1) for pt in self]))
        return values,

    def _table_title_str(self):
        return 'Points'

    @property
    def __item_class__(self):
        return Point

    def sort(self, key=attrgetter('x', 'y', 'z'), reverse=False):
        super().sort(key=key, reverse=reverse)

    def __cast(self, other):
        return other.data if isinstance(other, UserList) else other

    def _is_valid_operand(self, other):
        return isinstance(other, (list, np.ndarray, self.__class__,
                                  self.__item_class__))

    def _is_compatible_shape(self, other):
        other = np.asarray(self.__cast(other))
        return np.asarray(self.data).shape[-1] == len(other) if \
            other.ndim == 1 else np.asarray(self.data).shape == other.shape

    def _check_operands(self, other):
        if not self._is_compatible_shape(other):
            self_shape = np.asarray(self.data).shape
            other_shape = np.asarray(self.__cast(other)).shape
            raise ValueError(operand_shape_error_msg.format(self_shape,
                             other_shape))

    def __add__(self, other):
        if not ((self._is_valid_operand(other) and
                 self._is_compatible_shape(other)) or np.isscalar(other)):
            return NotImplemented
        if np.isscalar(other) or isinstance(other, self.__item_class__):
            data = [pt.__add__(other) for pt in self]
        else:
            other = np.asarray(self.__cast(other))
            if other.ndim == 1:
                return self.__add__(self.__item_class__(other))
            else:
                data = [pt.__add__(other_pt) for pt, other_pt in
                        zip(self.data, other)]
        return self.__class__(data)

    def __radd__(self, other):
        if not ((self._is_valid_operand(other) and
                 self._is_compatible_shape(other)) or np.isscalar(other)):
            return NotImplemented
        return self.__add__(other)

    def __iadd__(self, other):
        if not ((self._is_valid_operand(other) and
                 self._is_compatible_shape(other)) or np.isscalar(other)):
            return NotImplemented
        [self.__setitem__(i, pt.__iadd__(other)) for i, pt in
         enumerate(self)]
        return self

    def __sub__(self, other):
        if not ((self._is_valid_operand(other) and
                 self._is_compatible_shape(other)) or np.isscalar(other)):
            return NotImplemented
        if np.isscalar(other) or isinstance(other, self.__item_class__):
            data = [pt.__sub__(other) for pt in self]
        else:
            other = np.asarray(self.__cast(other))
            if other.ndim == 1:
                return self.__sub__(self.__item_class__(other))
            else:
                data = [pt.__sub__(other_pt) for pt, other_pt in
                        zip(self.data, other)]
        return self.__class__(data)

    def __rsub__(self, other):
        if not ((self._is_valid_operand(other) and
                 self._is_compatible_shape(other)) or np.isscalar(other)):
            return NotImplemented
        return self.__sub__(other)

    def __isub__(self, other):
        if not ((self._is_valid_operand(other) and
                 self._is_compatible_shape(other)) or np.isscalar(other)):
            return NotImplemented
        [self.__setitem__(i, pt.__isub__(other)) for i, pt in
         enumerate(self)]
        return self

    @property
    def A(self):
        """Return array of vectors."""
        return self.asarray()

    @property
    def T(self):
        """Return transpose of :class:`Points` as an \
            :class:`~numpy:numpy.ndarray`."""
        return self.asarray().T

    @property
    def M(self):
        """Return :class:`Points` as a :class:`~numpy:numpy.matrix`."""
        return self.asmatrix()

    def asarray(self):
        """Return :class:`Points` as an :class:`~numpy:numpy.ndarray`."""
        return np.asarray(self.tolist())

    def asmatrix(self):
        """Return :class:`Points` as a :class:`~numpy:numpy.matrix`."""
        return np.asmatrix(self.tolist())

    @property
    def x(self):
        """Return :math:`x` coordinates of `Point` objects as array."""
        return np.asarray([point.x for point in self])

    @x.setter
    def x(self, values):
        self._check_operands(values)
        [setattr(point, 'x', value) for point, value in zip(self, values)]

    @property
    def y(self):
        """Return :math:`y` coordinates of `Point` objects as array."""
        return np.asarray([point.y for point in self])

    @y.setter
    def y(self, values):
        self._check_operands(values)
        [setattr(point, 'y', value) for point, value in zip(self, values)]

    @property
    def z(self):
        """Return :math:`z` coordinates of `Point` objects as array."""
        return np.asarray([point.z for point in self])

    @z.setter
    def z(self, values):
        self._check_operands(values)
        [setattr(point, 'z', value) for point, value in zip(self, values)]

    def filter(self, condition, invert=False):
        """Filter `Points` by `condition`.

        Parameters
        ----------
        condition : array_like, bool
            Boolean index array having same shape as the initial dimensions
            of the list of `Points` being indexed.
        invert : bool, optional
            If `True`, the boolean array `condition` is inverted element-wise.

        Returns
        -------
        filtered_points : `Points`
            If `invert` is `False`, return the elements where `condition`
            is `True`.

            If `invert` is `True`, return the elements where `~condition`
            (i.e., numpy.invert(condition)) is `True`.

        """
        if invert:
            condition = ~condition
        return self.__class__(points=np.asarray(self)[condition].tolist())

    def rezero(self, epsilon=1.0e-10):
        """Set really really small coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        [point.rezero(epsilon=epsilon) for point in self]

    def rotate(self, angle=None, axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None,
               degrees=False, transform_matrix=None, verbose=False, **kwargs):
        """Rotate `Point`\ s coordinates.

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
        [point.rotate(transform_matrix=transform_matrix) for point in self]

    def translate(self, t):
        """Translate `Point`\ s by :class:`Vector` `t`.

        Parameters
        ----------
        v : :class:`Vector`

        """
        [point.translate(t) for point in self]

    def tolist(self):
        """Return `Points` as :class:`~python:list`"""
        # return np.asarray([pt.tolist() for pt in self]).tolist()
        return [pt.tolist() for pt in self]

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(points=self.tolist())
