# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes with velocity attributes (:mod:`sknano.core.atoms._velocity_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._velocity_atoms

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

__all__ = ['VelocityAtom', 'VelocityAtoms']


@total_ordering
class VelocityAtom(Atom):
    """An `Atom` class with velocity component attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    vx, vy, vz : float, optional
        :math:`v_x, v_y, v_z` components of `VelocityAtom` velocity.

    """
    def __init__(self, *args, vx=None, vy=None, vz=None, **kwargs):

        super().__init__(*args, **kwargs)

        self._v = Vector([vx, vy, vz])
        self.fmtstr = super().fmtstr + \
            ", vx={vx:.6f}, vy={vy:.6f}, vz={vz:.6f}"

    def __eq__(self, other):
        return self.v == other.v and super().__eq__(other)

    def __lt__(self, other):
        return (self.v < other.v and super().__le__(other)) or \
            (self.v <= other.v and super().__lt__(other))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['vx', 'vy', 'vz'])
        return attrs

    @property
    def vx(self):
        """:math:`x` component of `VelocityAtom` velocity vector"""
        return self._v.x

    @vx.setter
    def vx(self, value):
        """Set :math:`v_x`.

        Set :math:`v_x`, the :math:`x` component of `VelocityAtom` velocity
        vector.

        Parameters
        ----------
        value : float
            :math:`v_x` component of velocity

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._v.x = value

    @property
    def vy(self):
        """:math:`x` component of `VelocityAtom` velocity vector"""
        return self._v.y

    @vy.setter
    def vy(self, value):
        """Set :math:`v_y`.

        Set :math:`v_y`, the :math:`y` component of `VelocityAtom` velocity
        vector.

        Parameters
        ----------
        value : float
            :math:`v_y` component of velocity

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._v.y = value

    @property
    def vz(self):
        """:math:`z` component of `VelocityAtom` velocity vector"""
        return self._v.z

    @vz.setter
    def vz(self, value):
        """Set :math:`v_z`.

        Set :math:`v_z`, the :math:`z` component of `VelocityAtom` velocity
        vector.

        Parameters
        ----------
        value : float
            :math:`v_z` component of velocity

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._v.z = value

    @property
    def v(self):
        """:math:`v_x, v_y, v_z` array of velocity components.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`v_x`, :math:`v_y`, :math:`v_z`]
            velocity components of `VelocityAtom`.

        """
        return self._v

    @v.setter
    def v(self, value):
        """Set :math:`x, y, z` components of `VelocityAtom` velocity.

        Parameters
        ----------
        value : array_like
            3-element ndarray of [:math:`v_x`, :math:`v_y`, :math:`v_z`]
            velocity components of `VelocityAtom`.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._v[:] = Vector(value, nd=3)

    def rezero(self, epsilon=1.0e-10):
        """Re-zero velocity vector components.

        Set velocity vector components with absolute value less than
        `epsilon` to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        self.v.rezero(epsilon)
        super().rezero(epsilon)

    def rotate(self, **kwargs):
        """Rotate `Atom` velocity vector.

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
        self.v.rotate(**kwargs)
        super().rotate(**kwargs)

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(vx=self.vx, vy=self.vy, vz=self.vz))
        return super_dict


class VelocityAtoms(Atoms):
    """An `Atoms` class for `VelocityAtom`\ s.

    A container class for `VelocityAtom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `VelocityAtoms`}, optional
        if not `None`, then a list of `VelocityAtom` instance objects or an
        existing `VelocityAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return VelocityAtom

    def sort(self, key=attrgetter('v'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def velocities(self):
        """:class:`~numpy:numpy.ndarray` of `VelocityAtom` velocities."""
        return np.asarray([atom.v for atom in self])

    @property
    def v(self):
        """Alias for :attr:`~VelocityAtoms.velocities`."""
        return self.velocities

    @property
    def vx(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`v_x` \
            components."""
        return self.v[:, 0]

    @property
    def vy(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`v_y` \
            components."""
        return self.v[:, 1]

    @property
    def vz(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`v_z` \
            components."""
        return self.v[:, 2]
