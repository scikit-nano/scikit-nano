# -*- coding: utf-8 -*-
"""
==============================================================================
Custom NumPy Point class (:mod:`sknano.core.math._point`)
==============================================================================

.. currentmodule:: sknano.core.math._point

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numbers
import numpy as np

from ._transforms import rotate

__all__ = ['Point']


class Point(np.ndarray):
    """Abstract object representation for a point in :math:`R^n`.

    Parameters
    ----------
    p : array_like, optional
        :math:`x, y` coordinates of point in :math:`R^2` space.
        :math:`x, y, z` coordinates of point in :math:`R^3` space.
    nd : {None, int}, optional
    dtype : data-type, optional
    copy : bool, optional

    """
    __array_priority__ = 10.0

    def __new__(cls, p=None, nd=None, dtype=None, copy=True):

        if isinstance(p, Point):
            if nd is not None and isinstance(nd, numbers.Number) and \
                    len(p) < int(nd):
                p = np.append(p, np.zeros(int(nd) - len(p)))

            if dtype is None:
                intype = p.dtype
            else:
                intype = np.dtype(dtype)

            pt = p.view(cls)
            if intype != p.dtype:
                return pt.astype(intype)

            if copy:
                return pt.copy()
            else:
                return pt

        dtype = np.dtype(dtype)

        if isinstance(p, (tuple, list, np.ndarray)):
            try:
                for i, coord in enumerate(p[:]):
                    if coord is None:
                        p[i] = 0.0
            except TypeError:
                p = np.zeros(len(p), dtype=dtype)
            else:
                p = np.asarray(p, dtype=dtype)
            if nd is not None and isinstance(nd, numbers.Number) and \
                    len(p) < int(nd):
                p = np.append(p, np.zeros(int(nd) - len(p)))
            nd = len(p)
        else:
            if nd is None or not isinstance(nd, numbers.Number):
                nd = 3
            else:
                nd = int(nd)
            p = np.zeros(nd, dtype=dtype)

        arr = np.array(p, dtype=dtype, copy=copy).view(cls)
        pt = np.ndarray.__new__(cls, arr.shape, arr.dtype, buffer=arr)

        pt.nd = len(pt)

        return pt

    def __array_finalize__(self, pt):
        if pt is None:
            return None

        self.nd = len(pt)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return "Point({!r})".format(self.tolist())

    def tolist(self):
        """List of `Point` coordinates formatted for *pretty* output."""
        return np.around(self.__array__(), decimals=10).tolist()

    def __getattr__(self, name):
        try:
            nd = len(self)
            if nd == 2 and name in ('x', 'y'):
                if name == 'x':
                    return self[0]
                else:
                    return self[1]
            elif nd == 3 and name in ('x', 'y', 'z'):
                if name == 'x':
                    return self[0]
                elif name == 'y':
                    return self[1]
                else:
                    return self[2]
        except TypeError:
            pass
        return super().__getattribute__(name)

    def __setattr__(self, name, value):
        nd = getattr(self, 'nd', None)
        if nd is not None and nd == 2 and name in ('x', 'y'):
            if name == 'x':
                self[0] = value
            else:
                self[1] = value
        elif nd is not None and nd == 3 and name in ('x', 'y', 'z'):
            if name == 'x':
                self[0] = value
            elif name == 'y':
                self[1] = value
            else:
                self[2] = value
        else:
            super().__setattr__(name, value)

    def __eq__(self, other):
        if not isinstance(other, Point):
            other = Point(other)
        return self is other or np.allclose(self.__array__(),
                                            other.__array__())

    def __lt__(self, other):
        if not isinstance(other, Point):
            other = Point(other)
        # origin = Point(nd=self.nd)
        # return self.euclidean_distance(origin) < \
        #     other.euclidean_distance(origin)
        return np.all(np.less(self.__array__(), other.__array__())) or \
            (np.any(np.less(self.__array__(), other.__array__())) and
             (self.x < other.x) or (self.x <= other.x and self.y < other.y) or
             (self.x <= other.x and self.y <= other.y and self.z < other.z))

    def __le__(self, other):
        return self < other or self == other

    def __gt__(self, other):
        return not (self < other or self == other)

    def __ge__(self, other):
        return not (self < other)

    def __ne__(self, other):
        return not self == other

    @property
    def column_matrix(self):
        """Return column matrix representation of `Point` coordinates."""
        return np.matrix(self.__array__().reshape(self.shape[0], 1))

    @property
    def row_matrix(self):
        """Return row matrix representation of `Point` coordinates."""
        return np.matrix(self.__array__())

    def euclidean_distance(self, pt):
        """Compute the euclidean distance between `pt` and `self`."""
        return np.sqrt(((self.__array__() - pt.__array__()) ** 2).sum())

    def rezero_coords(self, epsilon=1.0e-10):
        """Alias for :meth:`Point.rezero`."""
        self.rezero(epsilon=epsilon)

    def rezero(self, epsilon=1.0e-10):
        """Re-zero `Point` coordinates near zero.

        Set `Point` coordinates with absolute value less than `epsilon` to
        zero.

        Parameters
        ----------
        epsilon : float, optional
            Smallest allowed absolute value of any :math:`x,y,z` coordinate.

        """
        self[np.where(np.abs(self.__array__()) <= epsilon)] = 0.0

    def rotate(self, angle=None, axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None, degrees=False,
               transform_matrix=None, verbose=False, **kwargs):
        """Rotate `Point` coordinates.

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
        core.math.rotate

        """
        self[:] = rotate(self, angle=angle, axis=axis,
                         anchor_point=anchor_point, rot_point=rot_point,
                         from_vector=from_vector, to_vector=to_vector,
                         transform_matrix=transform_matrix, degrees=degrees,
                         verbose=verbose, **kwargs)

    def translate(self, t):
        """Translate `Point` coordinates by :class:`~sknano.core.math.Vector` \
            `t`.

        Parameters
        ----------
        t : :class:`~sknano.core.math.Vector`

        See Also
        --------
        core.math.translate

        """
        self += t
