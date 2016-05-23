# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes with a dipole moment (:mod:`sknano.core.atoms.dipole_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms.dipole_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter
import numbers

# import numpy as np

from sknano.core.math import Vector, Vectors
from .atoms import Atom, Atoms

__all__ = ['DipoleAtom', 'DipoleAtoms']


class DipoleAtom(Atom):
    """An `Atom` sub-class with electric dipole moment attributes.

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

    # def __lt__(self, other):
    #     return (self.p < other.p and super().__le__(other)) or \
    #         (self.p <= other.p and super().__lt__(other))

    @property
    def __atoms_class__(self):
        return DipoleAtoms

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.p == other.p and super().__eq__(other)

    def __le__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.p > other.p or not super().__le__(other):
            return False
        return True

    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.p >= other.p or not super().__lt__(other):
            return False
        return True

    def __ge__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.p < other.p or not super().__ge__(other):
            return False
        return True

    def __gt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.p <= other.p or not super().__gt__(other):
            return False
        return True

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
        """Return :class:`~python:dict` of constructor parameters."""
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
    def p(self):
        """:class:`Vectors` of :attr:`DipoleAtom.p` :class:`Vector`\ s."""
        return Vectors([atom.p for atom in self])

    @property
    def dipole_moments(self):
        """An alias for :attr:`DipoleAtoms.p"""
        return self.p

    @property
    def px(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`DipoleAtom.px` \
            components."""
        return self.p.x

    @property
    def py(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`DipoleAtom.py` \
            components."""
        return self.p.y

    @property
    def pz(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`DipoleAtom.pz` \
            components."""
        return self.p.z

    @property
    def P(self):
        """Return the total net dipole moment vector of `DipoleAtoms`."""
        return Vector(self.dipole_moments.sum(axis=1))
