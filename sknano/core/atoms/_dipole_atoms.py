# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes with a dipole moment (:mod:`sknano.core.atoms._dipole_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._dipole_atoms

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

__all__ = ['DipoleAtom', 'DipoleAtoms']


@total_ordering
class DipoleAtom(Atom):
    """An `Atom` class with electric dipole moment attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    px, py, pz : float, optional
        :math:`p_x, p_y, p_z` components of `DipoleAtom` dipole moment.

    """
    def __init__(self, *args, px=None, py=None, pz=None, **kwargs):

        super().__init__(*args, **kwargs)

        self._p = Vector([px, py, pz])
        self.fmtstr = super().fmtstr + \
            ", px={px:.6f}, py={py:.6f}, pz={pz:.6f}"

    def __eq__(self, other):
        return self.p == other.p and super().__eq__(other)

    def __lt__(self, other):
        return (self.p < other.p and super().__le__(other)) or \
            (self.p <= other.p and super().__lt__(other))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['px', 'py', 'pz'])
        return attrs

    @property
    def px(self):
        """:math:`x` component of `DipoleAtom` dipole moment vector"""
        return self._p.x

    @px.setter
    def px(self, value):
        """Set :math:`p_x`.

        Set :math:`p_x`, the :math:`x` component of `DipoleAtom` dipole moment
        vector.

        Parameters
        ----------
        value : float
            :math:`p_x` component of dipole moment

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._p.x = value

    @property
    def py(self):
        """:math:`x` component of `DipoleAtom` dipole moment vector"""
        return self._p.y

    @py.setter
    def py(self, value):
        """Set :math:`p_y`.

        Set :math:`p_y`, the :math:`y` component of `DipoleAtom` dipole moment
        vector.

        Parameters
        ----------
        value : float
            :math:`p_y` component of dipole moment

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._p.y = value

    @property
    def pz(self):
        """:math:`z` component of `DipoleAtom` dipole moment vector"""
        return self._p.z

    @pz.setter
    def pz(self, value):
        """Set :math:`p_z`.

        Set :math:`p_z`, the :math:`z` component of `DipoleAtom` dipole moment
        vector.

        Parameters
        ----------
        value : float
            :math:`p_z` component of dipole moment

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._p.z = value

    @property
    def p(self):
        """Dipole moment :math:`\\mathbf{p}=q\\mathbf{d}`."""
        return self._p

    @p.setter
    def p(self, value):
        """Set `DipoleAtom` dipole moment :math:`\\mathbf{p}=q\\mathbf{d}`.

        Parameters
        ----------
        value : array_like
            electric dipole moment vector

        """
        self._p[:] = Vector(value, nd=3)

    def rezero(self, epsilon=1.0e-10):
        """Re-zero dipole moment vector components.

        Set dipole moment vector components with absolute value less than
        `epsilon` to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        self.p.rezero(epsilon)
        super().rezero(epsilon)

    def rotate(self, **kwargs):
        """Rotate `Atom` dipole moment vector.

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
        self.p.rotate(**kwargs)
        super().rotate(**kwargs)

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(px=self.px, py=self.py, pz=self.pz))
        return super_dict


class DipoleAtoms(Atoms):
    """An `Atoms` sub-class for `DipoleAtom`\ s.

    A container class for `DipoleAtom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `DipoleAtoms`}, optional
        if not `None`, then a list of `DipoleAtom` instance objects or an
        existing `DipoleAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return DipoleAtom

    def sort(self, key=attrgetter('p'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def dipole_moments(self):
        """Return array of `DipoleAtom` dipole moments."""
        return np.asarray([atom.p for atom in self])

    @property
    def P(self):
        """Return the total net dipole moment vector of `DipoleAtoms`."""
        return Vector(self.dipole_moments.sum(axis=1))
