# -*- coding: utf-8 -*-
"""
==============================================================================
Vectors class (:mod:`sknano.core.math._vectors`)
==============================================================================

.. currentmodule:: sknano.core.math._vectors

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from sknano.core import UserList
from ._transforms import transformation_matrix
# from sknano.core.geometric_regions import Cuboid  # , Rectangle

__all__ = ['Vectors']

operand_shape_error_msg = \
    "operands could not be broadcast together with shapes {}, {}"


class Vectors(UserList):
    """Container class for collection of `Vector` objects.

    Parameters
    ----------
    vectors : {None, sequence, `Vectors`}, optional
        if not `None`, then a list of `Vector` instance objects or an
        existing `Vectors` instance object.

    """

    def __init__(self, vectors=None):
        # if vectors is not None:
        #     vectors = np.asarray(vectors).tolist()
        super().__init__(initlist=vectors)
        self.fmtstr = "{vectors!r}"

    def __repr__(self):
        return str(np.asarray([vec.tolist() for vec in self]))

    @property
    def __item_class__(self):
        return Vector

    def sort(self, key=attrgetter('p0', 'length'), reverse=False):
        super().sort(key=key, reverse=reverse)

    def __cast(self, other):
        return other.data if isinstance(other, UserList) else other

    def _is_valid_operand(self, other):
        return isinstance(other, (tuple, list, np.ndarray, self.__class__,
                                  self.__item_class__))

    def _is_compatible_shape(self, other):
        other = np.asarray(self.__cast(other))
        return np.asarray(self.data).shape[-1] == len(other) if \
            other.ndim == 1 else np.asarray(self.data).shape == other.shape

    def _validate_operand(self, other):
        valid = self._is_valid_operand(other)
        compatible = self._is_compatible_shape(other)
        if not valid:
            raise TypeError('Invalid operand type {!r}.'.format(type(other)))
        if not compatible:
            self_shape = np.asarray(self.data).shape
            other_shape = np.asarray(self.__cast(other)).shape
            raise ValueError(operand_shape_error_msg.format(self_shape,
                             other_shape))

    def __add__(self, other):
        if not ((self._is_valid_operand(other) and
                 self._is_compatible_shape(other)) or np.isscalar(other)):
            return NotImplemented
        if np.isscalar(other) or isinstance(other, self.__item_class__):
            data = [vec.__add__(other) for vec in self]
        else:
            other = np.asarray(self.__cast(other))
            if other.ndim == 1:
                return self.__add__(self.__item_class__(other))
            else:
                data = [vec.__add__(other_vec) for vec, other_vec in
                        zip(self.data, other)]
        return self.__class__(data)

    def __radd__(self, other):
        if not ((self._is_valid_operand(other) and
                 self._is_compatible_shape(other)) or np.isscalar(other)):
            return NotImplemented
        return self.__add__(self.__item_class__(other))

    def __iadd__(self, other):
        if not ((self._is_valid_operand(other) and
                 self._is_compatible_shape(other)) or np.isscalar(other)):
            return NotImplemented
        [self.__setitem__(i, vec.__iadd__(other)) for i, vec in
         enumerate(self)]
        return self

    def __sub__(self, other):
        if not ((self._is_valid_operand(other) and
                 self._is_compatible_shape(other)) or np.isscalar(other)):
            return NotImplemented
        if np.isscalar(other) or isinstance(other, self.__item_class__):
            data = [vec.__sub__(other) for vec in self]
        else:
            other = np.asarray(self.__cast(other))
            if other.ndim == 1:
                return self.__sub__(self.__item_class__(other))
            else:
                data = [vec.__sub__(other_vec) for vec, other_vec in
                        zip(self.data, other)]
        return self.__class__(data)

    def __rsub__(self, other):
        if not ((self._is_valid_operand(other) and
                 self._is_compatible_shape(other)) or np.isscalar(other)):
            return NotImplemented
        return self.__add__(-self.__item_class__(other))

    def __isub__(self, other):
        if not ((self._is_valid_operand(other) and
                 self._is_compatible_shape(other)) or np.isscalar(other)):
            return NotImplemented
        [self.__setitem__(i, vec.__isub__(other)) for i, vec in
         enumerate(self)]
        return self

    def __pow__(self, other):
        if not np.isscalar(other):
            return NotImplemented
        data = [vec.__pow__(other) for vec in self]
        return self.__class__(data)

    # def __and__(self, other):
    #     if not self._is_valid_operand(other):
    #         return NotImplemented
    #     return self.___class__([atom for atom in self.__cast(other)
    #                             if atom in self])

    # __rand__ = __and__

    # def __iand__(self, other):
    #     if not self._is_valid_operand(other):
    #         return NotImplemented
    #     self.data -= list(atom for atom in self.__cast(other)
    #                       if atom not in self)
    #     return self

    # def __or__(self, other):
    #     if not self._is_valid_operand(other):
    #         return NotImplemented
    #     # data = self.data + list(atom for atom in other if atom not in self)
    #     data = dedupe((atom for data in (self, self.__cast(other))
    #                    for atom in data), key=attrgetter('id'))
    #     return self.__class__(data)

    # __ror__ = __or__

    # def __ior__(self, other):
    #     if not self._is_valid_operand(other):
    #         return NotImplemented
    #     self.data += \
    #         list(atom for atom in self.__cast(other) if atom not in self)
    #     return self

    # def __xor__(self, other):
    #     if not self._is_valid_operand(other):
    #         return NotImplemented
    #     other = self.__cast(other)
    #     return (self - other) | (other - self)

    # __rxor__ = __xor__

    # def __ixor__(self, other):
    #     if not self._is_valid_operand(other):
    #         return NotImplemented
    #     if other is self:
    #         self.clear()
    #     else:
    #         other = self.__cast(other)
    #         for atom in other:
    #             if atom in self:
    #                 self.remove(atom)
    #             else:
    #                 self.append(atom)
    #     return self

    def __neg__(self):
        vecs = self.__class__([-vec for vec in self])
        vecs.rezero()
        return vecs

    @property
    def A(self):
        """Return array of vectors."""
        return self.asarray()

    @property
    def T(self):
        """Return transpose of :class:`Vectors` as an \
            :class:`~numpy:numpy.ndarray`."""
        return self.asarray().T

    @property
    def M(self):
        """Return :class:`Vectors` as a :class:`~numpy:numpy.matrix`."""
        return self.asmatrix()

    def asarray(self):
        """Return :class:`Vectors` as an :class:`~numpy:numpy.ndarray`."""
        return np.asarray([vec.tolist() for vec in self])

    def asmatrix(self):
        """Return :class:`Vectors` as a :class:`~numpy:numpy.matrix`."""
        return np.asmatrix([vec.tolist() for vec in self])

    @property
    def x(self):
        """Return :math:`x` coordinates of `Vector` objects as array."""
        return np.asarray([vec.x for vec in self])

    @x.setter
    def x(self, values):
        self._validate_operand(values)
        [setattr(vec, 'x', val) for vec, val in zip(self, values)]

    @property
    def y(self):
        """Return :math:`y` coordinates of `Vector` objects as array."""
        return np.asarray([vec.y for vec in self])

    @y.setter
    def y(self, values):
        self._validate_operand(values)
        [setattr(vec, 'y', val) for vec, val in zip(self, values)]

    @property
    def z(self):
        """Return :math:`z` coordinates of `Vector` objects as array."""
        return np.asarray([vec.z for vec in self])

    @z.setter
    def z(self, values):
        self._validate_operand(values)
        [setattr(vec, 'z', val) for vec, val in zip(self, values)]

    def angle(self, other):
        """Angles between each :class:`Vector` with `other`."""
        self._validate_operand(other)
        if isinstance(other, self.__item_class__):
            return np.asarray([vec.angle(other) for vec in self])
        else:
            return np.asarray([vec.angle(other_vec) for vec, other_vec
                               in zip(self, self.__cast(other))])

    def cross(self, other):
        """Cross product of :class:`Vector`\ s with `other`."""
        self._validate_operand(other)
        if isinstance(other, self.__item_class__):
            return self.__class__([vec.cross(other) for vec in self])
        else:
            return self.__class__([vec.cross(other_vec) for vec, other_vec
                                   in zip(self, self.__cast(other))])

    def dot(self, other, out=None):
        """Dot product of :class:`Vector`\ s with `other`."""
        self._validate_operand(other)
        if isinstance(other, self.__item_class__):
            return np.asarray([vec.dot(other) for vec in self])
        else:
            return np.asarray([vec.dot(other_vec) for vec, other_vec
                               in zip(self, self.__cast(other))])

    def projection(self, other):
        """Vector projection of each :class:`Vector` onto `other`."""
        self._validate_operand(other)
        if isinstance(other, self.__item_class__):
            return self.__class__([vec.projection(other) for vec in self])
        else:
            return self.__class__([vec.projection(other_vec) for vec, other_vec
                                   in zip(self, self.__cast(other))])

    def rejection(self, v):
        """Vector rejection onto `other`."""
        u = self
        return u - self.projection(v)

    def normalize(self):
        """Normalize each :class:`Vector`."""
        [vec.normalize() for vec in self]

    @property
    def norms(self):
        """Return `Vector` :attr:`Vector.norm`\ s as array."""
        return np.asarray([v.norm for v in self])

    @property
    def lengths(self):
        """Alias for :attr:`~Vectors.norms`."""
        return self.norms

    @property
    def magnitudes(self):
        """Alias for :attr:`~Vectors.norms`."""
        return self.norms

    def filter(self, condition, invert=False):
        """Filter `Vectors` by `condition`.

        Parameters
        ----------
        condition : array_like, bool
            Boolean index array having same shape as the initial dimensions
            of the list of `Vectors` being indexed.
        invert : bool, optional
            If `True`, the boolean array `condition` is inverted element-wise.

        Returns
        -------
        filtered_vectors : `Vectors`
            If `invert` is `False`, return the elements where `condition`
            is `True`.

            If `invert` is `True`, return the elements where `~condition`
            (i.e., numpy.invert(condition)) is `True`.

        """
        if invert:
            condition = ~condition
        return self.__class__(vectors=np.asarray(self)[condition].tolist())

    def rezero(self, epsilon=1.0e-10):
        """Set really really small coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        [vector.rezero(epsilon=epsilon) for vector in self]

    def rotate(self, angle=None, axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None,
               degrees=False, transform_matrix=None,
               fix_anchor_points=False, verbose=False, **kwargs):
        """Rotate `Vector`\ s coordinates.

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
        if 'fix_anchor_point' in kwargs:
            fix_anchor_points = kwargs['fix_anchor_point']
            del kwargs['fix_anchor_point']

        if transform_matrix is None:
            transform_matrix = \
                transformation_matrix(angle=angle, axis=axis,
                                      anchor_point=anchor_point,
                                      rot_point=rot_point,
                                      from_vector=from_vector,
                                      to_vector=to_vector, degrees=degrees,
                                      verbose=verbose, **kwargs)
        [vector.rotate(fix_anchor_point=fix_anchor_points,
                       transform_matrix=transform_matrix) for vector in self]

    def scale(self):
        return NotImplemented

    def tolist(self):
        """Return `Vectors` as :class:`~python:list`"""
        return np.asarray([vec.tolist() for vec in self]).tolist()

    def translate(self, t, fix_anchor_points=False, **kwargs):
        """Translate `Vector`\ s by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_points : bool, optional

        """
        [vector.translate(t, fix_anchor_point=fix_anchor_points)
         for vector in self]

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(vectors=[vec.tolist() for vec in self])

from ._vector import Vector
