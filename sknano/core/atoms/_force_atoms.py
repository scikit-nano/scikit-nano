# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes with force attributes (:mod:`sknano.core.atoms._force_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._force_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from functools import total_ordering
from operator import attrgetter
import numbers
import numpy as np

from sknano.core.math import Vector
from ._atoms import Atom, Atoms

__all__ = ['ForceAtom', 'ForceAtoms']


@total_ordering
class ForceAtom(Atom):
    """An `Atom` class with force attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    fx, fy, fz : float, optional
        :math:`f_x, f_y, f_z` components of `ForceAtom` velocity.

    """
    def __init__(self, *args, fx=None, fy=None, fz=None, **kwargs):

        super().__init__(*args, **kwargs)

        self._f = Vector([fx, fy, fz])
        self.fmtstr = super().fmtstr + \
            ", fx={fx:.6f}, fy={fy:.6f}, fz={fz:.6f}"

    def __eq__(self, other):
        return self.f == other.f and super().__eq__(other)

    def __lt__(self, other):
        return (self.f < other.f and super().__le__(other)) or \
            (self.f <= other.f and super().__lt__(other))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['fx', 'fy', 'fz'])
        return attrs

    @property
    def fx(self):
        """:math:`x` component of `ForceAtom` force vector"""
        return self._f.x

    @fx.setter
    def fx(self, value):
        """Set :math:`f_x`.

        Set :math:`f_x`, the :math:`x` component of `ForceAtom` force vector.

        Parameters
        ----------
        value : float
            :math:`f_x` component of force vector.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._f.x = value

    @property
    def fy(self):
        """:math:`x` component of `ForceAtom` force vector"""
        return self._f.y

    @fy.setter
    def fy(self, value):
        """Set :math:`f_y`.

        Set :math:`f_y`, the :math:`y` component of `ForceAtom` force vector.

        Parameters
        ----------
        value : float
            :math:`f_y` component of force vector.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._f.y = value

    @property
    def fz(self):
        """:math:`z` component of `ForceAtom` force vector"""
        return self._f.z

    @fz.setter
    def fz(self, value):
        """Set :math:`f_z`.

        Set :math:`f_z`, the :math:`z` component of `ForceAtom` force vector.

        Parameters
        ----------
        value : float
            :math:`f_z` component of force vector.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._f.z = value

    @property
    def f(self):
        """:math:`f_x, f_y, f_z` array of force vector components.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`f_x`, :math:`f_y`, :math:`f_z`]
            force vector components of `ForceAtom`.

        """
        return self._f

    @f.setter
    def f(self, value):
        """Set :math:`x, y, z` components of `ForceAtom` force vector.

        Parameters
        ----------
        value : array_like
            3-element ndarray of [:math:`f_x`, :math:`f_y`, :math:`f_z`]
            force components of `ForceAtom`.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._f[:] = Vector(value, nd=3)

    def rezero(self, epsilon=1.0e-10):
        """Re-zero position vector components.

        Set position vector components with absolute value less than
        `epsilon` to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        self.f.rezero(epsilon)
        super().rezero(epsilon)

    def rotate(self, **kwargs):
        """Rotate `Atom` position vector.

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
        self.f.rotate(**kwargs)
        super().rotate(**kwargs)

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(fx=self.fx, fy=self.fy, fz=self.fz))
        return super_dict


class ForceAtoms(Atoms):
    """An `Atoms` sub-class for `ForceAtom`\ s.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.ForceAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `ForceAtoms`}, optional
        if not `None`, then a list of `ForceAtom` instance objects or an
        existing `ForceAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return ForceAtom

    def sort(self, key=attrgetter('f'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def forces(self):
        """:class:`~numpy:numpy.ndarray` of `ForceAtom` forces."""
        return np.asarray([atom.f for atom in self])

    @property
    def f(self):
        """Alias for :attr:`~ForceAtoms.forces`."""
        return self.forces

    @property
    def fx(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`f_x` \
            components."""
        return self.f[:, 0]

    @property
    def fy(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`f_y` \
            components."""
        return self.f[:, 1]

    @property
    def fz(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`f_z` \
            components."""
        return self.f[:, 2]
