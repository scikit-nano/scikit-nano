# -*- coding: utf-8 -*-
"""
==================================================================
Custom NumPy Vector class (:mod:`sknano.core.math._vector`)
==================================================================

.. currentmodule:: sknano.core.math._vector

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import copy
import numbers
import warnings

import numpy as np
np.seterr(all='warn')

from sknano.core import Singleton
from ._point import Point
from ._transforms import rotate, transformation_matrix

__all__ = ['Vector', 'angle', 'cross', 'dot', 'scalar_triple_product',
           'vector_triple_product', 'scalar_projection', 'vector_projection',
           'vector_rejection', 'projection', 'rejection',
           'e1', 'e2', 'e3', 'xhat', 'yhat', 'zhat',
           '_check_vector_compatibility', 'NullVector']


def _check_vector_compatibility(v1, v2):
    if len(v1) != len(v2):
        raise ValueError("{!r} and {!r} must have same number "
                         "of components".format(v1, v2))


class Vector(np.ndarray):
    """Abstract object representation of a vector in :math:`R^n`

    Parameters
    ----------
    v : array_like, optional
        components of vector
    nd : {None, int}, optional
    p0 : array_like, optional
        Origin `Point` of vector in :math:`R^n` space.
    p : array_like, optional
        Terminating `Point` of vector.
        :math:`x, y` coordinates of point in :math:`R^2` space.
        :math:`x, y, z` coordinates of point in :math:`R^3` space.
    dtype : data-type, optional
    copy : bool, optional

    Notes
    -----
    .. todo::

       add new methods for coordinate transformations

    Examples
    --------
    I began writing my own `Point` and `Vector` classes while trying to
    teach myself about subclassing :class:`~numpy:numpy.ndarray`.
    I still don't completely understand the machinery of it, but I ended
    up with some code that's been useful for handling math operations involving
    points and vectors. More testing is in order...

    Here are some examples of their use.

    """
    __array_priority__ = 15.0

    def __new__(cls, v=None, nd=None, p0=None, p=None, dtype=None, copy=True):

        if isinstance(v, Vector):
            if nd is not None and isinstance(nd, numbers.Number) and \
                    len(v) < int(nd):
                v = np.append(v, np.zeros(int(nd) - len(v)))

            if dtype is None:
                intype = v.dtype
            else:
                intype = np.dtype(dtype)

            vec = v.view(cls)

            if p0 is not None:
                vec = Vector(np.asarray(vec), nd=nd,
                             p0=Point(p0, nd=nd, dtype=dtype, copy=copy))

            if intype != v.dtype:
                return vec.astype(intype)

            if copy:
                return vec.copy()
            else:
                return vec

        dtype = np.dtype(dtype)

        if isinstance(v, (tuple, list, np.ndarray)):
            try:
                for i, coord in enumerate(v[:]):
                    if coord is None:
                        v[i] = 0.0
            except TypeError:
                v = np.zeros(len(v), dtype=dtype)
            else:
                v = np.asarray(v, dtype=dtype)

            if nd is not None and isinstance(nd, numbers.Number) and \
                    len(v) < int(nd):
                v = np.append(v, np.zeros(int(nd) - len(v)))
            nd = len(v)

            if p0 is None:
                p0 = Point(nd=nd, dtype=dtype)
            else:
                p0 = Point(p0, nd=nd, dtype=dtype, copy=copy)
            p = p0 + v
        else:
            if p is None and p0 is None and \
                    (nd is None or not isinstance(nd, numbers.Number)):
                nd = 3
            if p is None:
                p = Point(nd=nd, dtype=dtype)
            else:
                p = Point(p, nd=nd, dtype=dtype, copy=copy)

            if p0 is None:
                p0 = Point(nd=nd, dtype=dtype)
            else:
                p0 = Point(p0, nd=nd, dtype=dtype, copy=copy)

            v = p - p0

        arr = np.array(v, dtype=dtype, copy=copy).view(cls)
        vec = np.ndarray.__new__(cls, arr.shape, arr.dtype, buffer=arr)

        vec.nd = len(vec)
        vec._p = p
        vec._p0 = p0

        return vec

    def __array_finalize__(self, obj):
        if obj is None:
            return None

        self.nd = len(obj)
        self._p0 = getattr(obj, 'p0', None)
        self._p = getattr(obj, 'p', None)
        if self._p0 is not None and self._p is None:
            try:
                self._p = self._p0 + self.__array__()
            except TypeError:
                try:
                    self._p = self._p0 + np.asarray(obj)
                except TypeError:
                    pass

    # def __array_prepare__(self, obj, context=None):
    #     if self.__array_priority__ >= Vector.__array_priority__:
    #         res = obj if isinstance(obj, type(self)) \
    #             else obj.view(type(self))
    #     else:
    #         res = obj.view(Vector)

    #     if context is None:
    #         return res
    #     return super(Vector, self).__array_prepare__(obj, context)

    def __array_wrap__(self, obj, context=None):
        res = np.ndarray.__array_wrap__(self, obj, context)
        return self.__class__(res.__array__(), p0=self.p0)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        try:
            if np.allclose(self.p0, np.zeros_like(self.p0)):
                return "Vector({!r})".format(self.tolist())
            else:
                return "Vector({!r}, p0={!r}, p={!r})".format(
                    self.tolist(), self.p0.tolist(), self.p.tolist())
        except AttributeError:
            return "Vector({!r})".format(self.tolist())

    def tolist(self):
        """List of `Vector` coordinates with values formatted for *pretty* \
            output."""
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
        return np.ndarray.__getattribute__(self, name)

    def __setattr__(self, name, value):
        # nd = len(self)
        nd = getattr(self, 'nd', None)

        if nd is not None and nd == 2 and name in ('x', 'y'):
            if name == 'x':
                self[0] = value
                try:
                    self._p.x = self._p0.x + value
                except AttributeError:
                    pass
            else:
                self[1] = value
                try:
                    self._p.y = self._p0.y + value
                except AttributeError:
                    pass
        elif nd is not None and nd == 3 and name in ('x', 'y', 'z'):
            if name == 'x':
                self[0] = value
                try:
                    self._p.x = self._p0.x + value
                except AttributeError:
                    pass
            elif name == 'y':
                self[1] = value
                try:
                    self._p.y = self._p0.y + value
                except AttributeError:
                    pass
            else:
                self[2] = value
                try:
                    self._p.z = self._p0.z + value
                except AttributeError:
                    pass
        else:
            np.ndarray.__setattr__(self, name, value)
            # if name is not None and name.endswith('p'):
            #    try:
            #        self[:] = self._p.__array__() - self._p0.__array__()
            #    except (AttributeError, TypeError):
            #        pass

    def __getitem__(self, index):
        data = np.ndarray.__getitem__(np.ndarray.view(self, np.ndarray),
                                      index)
        p0 = np.ndarray.__getitem__(np.ndarray.view(
                                    np.ndarray.__getattribute__(self, 'p0'),
                                    np.ndarray), index)
        p = np.ndarray.__getitem__(np.ndarray.view(
                                   np.ndarray.__getattribute__(self, 'p'),
                                   np.ndarray), index)

        try:
            data = data.view(type(self))
            data._p0 = np.ndarray.view(p0, Point)
            data._p = np.ndarray.view(p, Point)
            data._nd = len(data)
        except (AttributeError, TypeError):
            pass
        return data

    def __setitem__(self, index, value):
        data = np.ndarray.view(self, np.ndarray)
        p0 = np.ndarray.view(np.ndarray.__getattribute__(self, 'p0'),
                             np.ndarray)
        p = np.ndarray.view(np.ndarray.__getattribute__(self, 'p'),
                            np.ndarray)

        np.ndarray.__setitem__(data, index, value)
        np.ndarray.__setitem__(p, index, np.ndarray.__getitem__(p0, index) +
                               np.ndarray.__getitem__(data, index))

        data = data.view(type(self))
        data._p0 = np.ndarray.view(p0, Point)
        data._p = np.ndarray.view(p, Point)

    def __eq__(self, other):
        if not isinstance(other, Vector):
            other = Vector(other)
        return self is other or (np.allclose(self.__array__(),
                                             other.__array__())
                                 and np.allclose(self.p0, other.p0)
                                 and np.allclose(self.p, other.p))

    def __lt__(self, other):
        if not isinstance(other, Vector):
            other = Vector(other)
        return self.norm < other.norm

    def __le__(self, other):
        return self < other or self == other

    def __gt__(self, other):
        return not (self < other or self == other)

    def __ge__(self, other):
        return not (self < other)

    def __ne__(self, other):
        return not (self == other)

    def __mul__(self, other):
        if np.isscalar(other):
            return self.__class__(self.__array__() * other, p0=self.p0)
        elif isinstance(other, Vector) and other.nd == self.nd:
            print("Computing *scalar product* of Vector's:\n"
                  "{!r}\n{!r}".format(self, other))
            return self.dot(other)
        # elif isinstance(other, np.matrix):
        #     res = self.row_matrix * other
        #     if len(self) == res.shape[1]:
        #         return self.__class__(res.A.flatten(), p0=self.p0)
        #     elif res.shape == (1, 1):
        #         return res.A.flatten()[0]
        return NotImplemented

    def __rmul__(self, other):
        if np.isscalar(other):
            return self.__class__(other * self.__array__(), p0=self.p0)
        # elif isinstance(other, np.matrix):
        #     res = other * self.column_matrix
        #     if len(self) == res.shape[0]:
        #         return self.__class__(res.A.flatten(), p0=self.p0)
        #     elif res.shape == (1, 1):
        #         return res.A.flatten()[0]
        return NotImplemented

    def __truediv__(self, other):
        if np.isscalar(other):
            return Vector(self.__array__() / other, p0=self.p0)
        return NotImplemented

    __div__ = __truediv__

    def __floordiv__(self, other):
        if np.isscalar(other):
            return Vector(self.__array__() // other, p0=self.p0)
        return NotImplemented

    def __pow__(self, other, *modulo):
        if isinstance(other, numbers.Number):
            return Vector(self.__array__() ** other, p0=self.p0)
        return NotImplemented

    def __iadd__(self, other):
        """Add other to self in-place."""
        super().__iadd__(other)
        self._update_p()
        return self

    def __isub__(self, other):
        """Subtract other from self in-place."""
        super().__isub__(other)
        self._update_p()
        return self

    def __imul__(self, other):
        """Multiply self by other in-place."""
        if np.isscalar(other):
            super().__imul__(other)
            self._update_p()
            return self
        return NotImplemented

    def __itruediv__(self, other):
        if np.isscalar(other):
            super().__itruediv__(other)
            self._update_p()
            return self
        return NotImplemented

    __idiv__ = __itruediv__

    def __ifloordiv__(self, other):
        if np.isscalar(other):
            super().__ifloordiv__(other)
            self._update_p()
            return self
        return NotImplemented

    def __ipow__(self, other):
        if np.isscalar(other):
            super().__ipow__(other)
            self._update_p()
            return self
        return NotImplemented

    def copy(self):
        """Return a copy of the `Vector`."""
        return self.__copy__()

    def __copy__(self):
        try:
            return self.__class__(self.__array__(), p0=self.p0.__array__())
        except AttributeError:
            return self.__class__(self.__array__())

    def __deepcopy__(self, memo):
        return self.__copy__()

    def _update_p(self):
        self._p[:] = self._p0[:] + self.__array__()

    @property
    def length(self):
        """Alias for :attr:`Vector.norm`."""
        return self.norm

    @property
    def magnitude(self):
        """Alias for :attr:`Vector.norm`."""
        return self.norm

    @property
    def mag(self):
        """Alias for :attr:`Vector.norm`."""
        return self.norm

    @property
    def norm(self):
        """Return the vector norm."""
        return np.sqrt((self.__array__() ** 2).sum())

    @property
    def unit_vector(self):
        """Return unit vector."""
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            return self / self.norm

    @property
    def p(self):
        """:class:`Point` of `Vector` *head*."""
        return self._p

    @p.setter
    def p(self, value=np.ndarray):
        """Set new terminating :class:`Point` of `Vector`."""
        self._p[:] = value
        self[:] = self._p - self._p0

    @property
    def p0(self):
        """:class:`Point` of `Vector` *tail*."""
        return self._p0

    @p0.setter
    def p0(self, value=np.ndarray):
        """Set new origin :class:`Point` of vector."""
        self._p0[:] = value
        self[:] = self._p - self._p0

    def _translate_p0(self, t, fix_components=False):
        if fix_components:
            self.translate(t)
        else:
            self.p0.translate(t)

    def _translate_p(self, t, fix_components=False):
        if fix_components:
            self.translate(t)
        else:
            self.p.translate(t)

    @property
    def column_matrix(self):
        """Return column matrix representation of `Vector` coordinates."""
        return np.matrix(self.__array__().reshape(self.shape[0], 1))

    @property
    def row_matrix(self):
        """Return row matrix representation of `Vector` coordinates."""
        return np.matrix(self.__array__())

    def angle(self, other):
        """Angle between two `Vector`\ s."""
        _check_vector_compatibility(self, other)
        return np.arccos(np.dot(self.__array__(), other.__array__()) /
                         (self.norm * other.norm))

    def cross(self, other):
        """Cross product of two `Vector`\ s."""
        _check_vector_compatibility(self, other)
        val = np.cross(self.__array__(), other.__array__())
        if val.shape == ():
            return val[()]
        return Vector(val, p0=self.p0)

    def dot(self, other, out=None):
        """Dot product of two `Vector`\ s."""
        _check_vector_compatibility(self, other)
        return self.__array__().dot(other.__array__())

    def normalize(self):
        """Normalize the `Vector` to a :attr:`~Vector.unit_vector`."""
        self[:] = self.unit_vector

    def projection(self, v):
        """Vector projection onto `Vector` `v`."""
        u = self
        return dot(u, v) / dot(v, v) * v

    def rejection(self, v):
        """Vector rejection onto `Vector` `v`."""
        u = self
        return u - self.projection(v)

    def rezero_components(self, epsilon=1.0e-10):
        """Alias for :meth:`Vector.rezero`"""
        self.rezero(epsilon=epsilon)

    def rezero(self, epsilon=1.0e-10):
        """Re-zero `Vector` coordinates near zero.

        Set `Vector` coordinates with absolute value less than `epsilon` to
        zero.

        Parameters
        ----------
        epsilon : float, optional
            Smallest allowed absolute value of any :math:`x,y,z` coordinate.

        """
        self[np.where(np.abs(self.__array__()) <= epsilon)] = 0.0

    def rotate(self, angle=None, axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None,
               degrees=False, transform_matrix=None,
               fix_anchor_point=False, verbose=False, **kwargs):
        """Rotate `Vector` coordinates.

        Parameters
        ----------
        angle : float
            Rotation angle in radians, unless `degrees` is `True`.
        axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        degrees : :class:`~python:bool`, optional
            If `True`, you are saying that the `angle` is in degrees.
        transform_matrix : :class:`~numpy:numpy.ndarray`
        fix_anchor_point : :class:`~python:bool`, optional
            If `True`, leave the *tail* of the vector (:attr:`Vector.p0`)
            fixed (default: `False`).

        See Also
        --------
        core.math.rotate

        """
        if transform_matrix is None:
            transform_matrix = \
                transformation_matrix(angle=angle, axis=axis,
                                      anchor_point=anchor_point,
                                      rot_point=rot_point,
                                      from_vector=from_vector,
                                      to_vector=to_vector, degrees=degrees,
                                      verbose=verbose, **kwargs)

        self.p = rotate(self.p, transform_matrix=transform_matrix)
        if not fix_anchor_point:
            self.p0 = rotate(self.p0, transform_matrix=transform_matrix)

    def scale(self):
        return NotImplemented

    def translate(self, t, fix_anchor_point=False):
        """Translate `Vector` coordinates by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_point : bool, optional

        See Also
        --------
        core.math.translate

        """
        self.p += t
        if not fix_anchor_point:
            self.p0 += t


def angle(u, v):
    """Compute the angle between two Cartesian `Vector`\ s.

    Parameters
    ----------
    u, v : `Vector`

    Returns
    -------
    :class:`~numpy:numpy.number`

    """
    return np.arccos(dot(u, v) / (u.norm * v.norm))


def cross(u, v, p0=None):
    """Vector cross product of two `Vector`\ s.

    Parameters
    ----------
    u, v : `Vector`
    p0 : `Point`, optional

    Returns
    -------
    :class:`~numpy:numpy.number` or :class:`Vector`

    """
    val = np.cross(np.asarray(u), np.asarray(v))
    if p0 is None:
        p0 = u.p0
    if val.shape == ():
        return val[()]
    else:
        return Vector(val, p0=p0)


def dot(u, v):
    """Dot product of two `Vector`\ s.

    Parameters
    ----------
    u, v : `Vector`

    Returns
    -------
    :class:`~numpy:numpy.number`

    """
    return np.dot(np.asarray(u), np.asarray(v))


def scalar_triple_product(u, v, w):
    """Compute scalar triple product of three `Vector`\ s.

    Parameters
    ----------
    u, v, w : `Vector`

    Returns
    -------
    :class:`~numpy:numpy.number`

    """
    return dot(u, cross(v, w))


def vector_triple_product(u, v, w):
    """Compute vector triple product of three `Vector`\ s.

    Parameters
    ----------
    u, v, w : `Vector`

    Returns
    -------
    :class:`Vector`

    """

    return cross(u, cross(v, w))


def scalar_projection(a, b):
    """Compute the scalar projection of :math:`\\mathbf{a}` onto \
        :math:`\\mathbf{b}`.

    Parameters
    ----------
    a, b : `Vector`

    Returns
    -------
    :class:`~numpy:numpy.number`

    """
    return dot(a, b) / b.norm


def vector_projection(a, b):
    """Compute the vector projection of :math:`\\mathbf{a}` onto \
        :math:`\\mathbf{b}`.

    Parameters
    ----------
    a, b : `Vector`

    Returns
    -------
    :class:`Vector`

    """
    return dot(a, b) / dot(b, b) * b

projection = vector_projection


def vector_rejection(a, b):
    """Compute the vector rejection of :math:`\\mathbf{a}` onto \
        :math:`\\mathbf{b}`.

    Parameters
    ----------
    a, b : `Vector`

    Returns
    -------
    :class:`Vector`

    """
    a1 = vector_projection(a, b)
    return a - a1

rejection = vector_rejection

e1 = xhat = Vector([1, 0, 0])
e2 = yhat = Vector([0, 1, 0])
e3 = zhat = Vector([0, 0, 1])


class NullVector(metaclass=Singleton):
    def __init__(self):
        self.__instance = Vector([0, 0, 0])

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__,
                               self.__instance.tolist())
