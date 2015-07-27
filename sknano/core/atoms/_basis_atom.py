# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class for crystal structure basis (:mod:`sknano.core.atoms._basis_atom`)
===============================================================================

An `Atom` sub-class for a crystal structure basis atom.

.. currentmodule:: sknano.core.atoms._basis_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from functools import total_ordering
import copy
import numbers
import numpy as np

from sknano.core.math import Vector
from ._xyz_atom import XYZAtom

__all__ = ['BasisAtom']


@total_ordering
class BasisAtom(XYZAtom):
    """An abstract object representation of a crystal structure basis atom.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.Crystal3DLattice`
    xs, ys, zs : float
    """

    def __init__(self, *args, lattice=None, xs=None, ys=None, zs=None,
                 **kwargs):

        super().__init__(*args, **kwargs)

        # if None in (xs, ys, zs) and lattice is not None:
        #     xs, ys, zs = lattice.cartesian_to_fractional

        self.lattice = lattice
        try:
            self.rs = Vector([xs, ys, zs])
            # print(self.rs)
            # print(self.r)
        except AttributeError:
            pass
        self.fmtstr = super().fmtstr + \
            ", xs={xs:.6f}, ys={ys:.6f}, zs={zs:.6f}, lattice={lattice!r}"

    def __eq__(self, other):
        return self.rs == other.rs and super().__eq__(other)

    def __lt__(self, other):
        """Test if `self` is *less than* `other`."""
        return (self.rs < other.rs and super().__le__(other)) or \
            (self.rs <= other.rs and super().__lt__(other))

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['lattice', 'xs', 'ys', 'zs'])
        return attrs

    @property
    def xs(self):
        """:math:`x`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`x`-coordinate in units of **Angstroms**.

        """
        return self.rs.x

    @xs.setter
    def xs(self, value):
        """Set `Atom` :math:`x`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`x`-coordinate in units of **Angstroms**.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        rs = self.rs
        rs.x = value
        self._update_cartesian_coordinate(rs)

    @property
    def ys(self):
        """:math:`y`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        return self.rs.y

    @ys.setter
    def ys(self, value):
        """Set `Atom` :math:`y`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        rs = self.rs
        rs.y = value
        self._update_cartesian_coordinate(rs)

    @property
    def zs(self):
        """:math:`z`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        return self.rs.z

    @zs.setter
    def zs(self, value):
        """Set `Atom` :math:`z`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        rs = self.rs
        rs.z = value
        self._update_cartesian_coordinate(rs)

    @property
    def rs(self):
        """:math:`x, y, z` components of `Atom` position vector.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates of `Atom`.

        """
        try:
            return Vector(self.lattice.cartesian_to_fractional(self.r))
        except AttributeError:
            return Vector()

    @rs.setter
    def rs(self, value):
        """Set :math:`x, y, z` components of `Atom` position vector.

        Parameters
        ----------
        value : array_like
            :math:`x, y, z` coordinates of `Atom` position vector relative to
            the origin.

        """
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('Expected an array_like object')
        rs = Vector(value, nd=3)
        self._update_cartesian_coordinate(rs)

    def _update_cartesian_coordinate(self, rs):
        self.r = Vector(self.lattice.fractional_to_cartesian(rs))

    @property
    def lattice(self):
        """:class:`~sknano.core.crystallography.Crystal3DLattice`."""
        return self._lattice

    @lattice.setter
    def lattice(self, value):
        self._lattice = copy.copy(value)

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
        self.lattice.rotate(**kwargs)
        super().rotate(**kwargs)

    def translate(self, t, fix_anchor_point=True):
        """Translate `Atom` position vector by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_point : bool, optional

        """
        # TODO compare timing benchmarks for translation of position vector.
        self.lattice.translate(t, fix_anchor_point=fix_anchor_point)
        super().translate(t, fix_anchor_point=fix_anchor_point)

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(lattice=self.lattice, xs=self.xs, ys=self.ys,
                               zs=self.zs))
        return super_dict
