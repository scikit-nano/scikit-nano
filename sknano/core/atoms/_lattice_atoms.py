# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes for crystal lattices (:mod:`sknano.core.atoms._lattice_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._lattice_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# from operator import attrgetter
import copy
import numbers

import numpy as np

from sknano.core.math import Vector, Vectors

from ._atoms import Atom, Atoms
from .mixins import PBCAtomsMixin

__all__ = ['LatticeAtom', 'LatticeAtoms']


class LatticeAtom(Atom):
    """An `Atom` sub-class with crystal lattice attributes.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.Crystal3DLattice`
    xs, ys, zs : :class:`~python:float`

    """

    def __init__(self, *args, lattice=None, xs=None, ys=None, zs=None,
                 **kwargs):

        super().__init__(*args, **kwargs)

        self.lattice = lattice

        if all([x is not None for x in (xs, ys, zs)]):
            self.rs = Vector([xs, ys, zs])

        self.fmtstr = super().fmtstr + \
            ", lattice={lattice!r}, xs={xs!r}, ys={ys!r}, zs={zs!r}"
        # ", lattice={lattice!r}, xs={xs:.6f}, ys={ys:.6f}, zs={zs:.6f}"

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.rs is None or other.rs is None:
            return super().__eq__(other)
        return self.rs == other.rs and super().__eq__(other)

    def __le__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.rs is None or other.rs is None:
            return super().__le__(other)
        if self.rs > other.rs or not super().__le__(other):
            return False
        return True

    # def __lt__(self, other):
    #     """Test if `self` is *less than* `other`."""
    #     return (self.rs < other.rs and super().__le__(other)) or \
    #         (self.rs <= other.rs and super().__lt__(other))

    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.rs is None or other.rs is None:
            return super().__lt__(other)
        # return ((self.rs < other.rs and self.__le__(other)) or
        #         (self.__le__(other) and super().__lt__(other)))
        if self.rs >= other.rs or not super().__lt__(other):
            return False
        return True

    def __ge__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.rs is None or other.rs is None:
            return super().__ge__(other)
        if self.rs < other.rs or not super().__ge__(other):
            return False
        return True

    def __gt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.rs is None or other.rs is None:
            return super().__gt__(other)
        if self.rs <= other.rs or not super().__gt__(other):
            return False
        return True

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['lattice', 'xs', 'ys', 'zs'])
        return attrs

    @property
    def xs(self):
        """Scaled :math:`x`-coordinate.

        Returns
        -------
        float
            Scaled :math:`x`-coordinate.

        """
        try:
            return self.rs.x
        except AttributeError:
            return None

    @xs.setter
    def xs(self, value):
        """Set `Atom` :math:`x`-coordinate.

        Parameters
        ----------
        value : float
            :math:`x`-coordinate..

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')

        try:
            rs = self.rs
            rs.x = value
            self._update_cartesian_coordinate(rs)
        except AttributeError:
            pass

    @property
    def ys(self):
        """Scaled :math:`y`-coordinate.

        Returns
        -------
        float
            Scaled :math:`y`-coordinate.

        """
        try:
            return self.rs.y
        except AttributeError:
            return None

    @ys.setter
    def ys(self, value):
        """Set `Atom` :math:`y`-coordinate.

        Parameters
        ----------
        value : float
            :math:`y`-coordinate.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        try:
            rs = self.rs
            rs.y = value
            self._update_cartesian_coordinate(rs)
        except AttributeError:
            pass

    @property
    def zs(self):
        """Scaled :math:`z`-coordinate.

        Returns
        -------
        float
            Scaled :math:`z`-coordinate.

        """
        try:
            return self.rs.z
        except AttributeError:
            return None

    @zs.setter
    def zs(self, value):
        """Set `Atom` :math:`z`-coordinate.

        Parameters
        ----------
        value : float
            Scaled :math:`z`-coordinate

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        try:
            rs = self.rs
            rs.z = value
            self._update_cartesian_coordinate(rs)
        except AttributeError:
            pass

    @property
    def rs(self):
        """Scaled :math:`x, y, z` components of `Atom` position vector.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x_s, y_s, z_s`] coordinates of `Atom`.

        """
        try:
            return Vector(self.lattice.cartesian_to_fractional(self.r))
        except AttributeError:
            return None

    @rs.setter
    def rs(self, value):
        """Set scaled :math:`x, y, z` components of `Atom` position vector.

        Parameters
        ----------
        value : array_like
            :math:`x, y, z` coordinates of `Atom` position vector relative to
            the origin.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        self._update_cartesian_coordinate(Vector(value, nd=3))

    def _update_cartesian_coordinate(self, rs):
        try:
            self.r = self.lattice.fractional_to_cartesian(rs)
        except AttributeError:
            pass

    @property
    def lattice(self):
        """:class:`~sknano.core.crystallography.Crystal3DLattice`."""
        return self._lattice

    @lattice.setter
    def lattice(self, value):
        self._lattice = copy.deepcopy(value)
        try:
            self.rs = self.lattice.cartesian_to_fractional(self.r)
        except AttributeError:
            pass

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        super_dict = super().todict()
        super_dict.update(dict(lattice=self.lattice, xs=self.xs, ys=self.ys,
                               zs=self.zs))
        return super_dict


class LatticeAtoms(PBCAtomsMixin, Atoms):
    """An `Atoms` sub-class for crystal structure lattice atoms.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.LatticeAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `LatticeAtoms`}, optional
        if not `None`, then a list of `LatticeAtom` instance objects or an
        existing `LatticeAtoms` instance object.

    """

    @property
    def __atom_class__(self):
        return LatticeAtom

    @property
    def rs(self):
        """:class:`Vectors` of :attr:`LatticeAtom.rs` :class:`Vector`\ s"""
        return Vectors([atom.rs for atom in self])

    @property
    def xs(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`LatticeAtom.xs` values"""
        return self.rs.x

    @property
    def ys(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`LatticeAtom.ys` values"""
        return self.rs.y

    @property
    def zs(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`LatticeAtom.zs` values"""
        return self.rs.z

    @property
    def lattice(self):
        """Return the :attr:`LatticeAtom.lattice` of the first atom in self."""
        try:
            return self[0].lattice
        except IndexError:
            return None

    @lattice.setter
    def lattice(self, value):
        [setattr(atom, 'lattice', value) for atom in self]

    @property
    def cell_matrix(self):
        """Return the :attr:`Crystal3DLattice.cell_matrix`."""
        try:
            return self.lattice.cell_matrix
        except AttributeError:
            return None

    @property
    def cell(self):
        """Alias for :attr:`LatticeAtoms.cell_matrix`."""
        return self.cell_matrix

    def wrap_coords(self, pbc=None):
        """Wrap coordinates into lattice."""
        try:
            [setattr(atom, 'r', self.lattice.wrap_cartesian_coordinate(
                     atom.r, pbc=pbc if pbc is not None else self.pbc))
             for atom in self]
        except AttributeError:
            pass
