# -*- coding: utf-8 -*-
"""
===============================================================================
Lattice molecule classes (:mod:`sknano.core.molecules._lattice_molecules`)
===============================================================================

.. currentmodule:: sknano.core.molecules._lattice_molecules

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# from operator import attrgetter
from functools import total_ordering
import copy
import numbers
import numpy as np

from sknano.core.math import Vector

from ._molecules import Molecule, Molecules

__all__ = ['LatticeMolecule', 'LatticeMolecules']


@total_ordering
class LatticeMolecule(Molecule):
    """Class representation of a crystal structure lattice molecule.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.Crystal3DLattice`
    xs, ys, zs : float
    """

    def __init__(self, *args, lattice=None, xs=None, ys=None, zs=None,
                 **kwargs):

        super().__init__(*args, **kwargs)

        self.lattice = lattice

        if all([x is not None for x in (xs, ys, zs)]):
            self.rs = Vector([xs, ys, zs])

        self.fmtstr = super().fmtstr + \
            ", lattice={lattice!r}, xs={xs!r}, ys={ys!r}, zs={zs!r}"

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
        try:
            return self.rs.x
        except AttributeError:
            return None

    @xs.setter
    def xs(self, value):
        """Set `Molecule` :math:`x`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`x`-coordinate in units of **Angstroms**.

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
        """:math:`y`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`y`-coordinate in units of **Angstroms**.

        """
        try:
            return self.rs.y
        except AttributeError:
            return None

    @ys.setter
    def ys(self, value):
        """Set `Molecule` :math:`y`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`y`-coordinate in units of **Angstroms**.

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
        """:math:`z`-coordinate in units of **Angstroms**.

        Returns
        -------
        float
            :math:`z`-coordinate in units of **Angstroms**.

        """
        try:
            return self.rs.z
        except AttributeError:
            return None

    @zs.setter
    def zs(self, value):
        """Set `Molecule` :math:`z`-coordinate in units of **Angstroms**.

        Parameters
        ----------
        value : float
            :math:`z`-coordinate in units of **Angstroms**.

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
        """:math:`x, y, z` components of `Molecule` position vector.

        Returns
        -------
        ndarray
            3-element ndarray of [:math:`x, y, z`] coordinates of `Molecule`.

        """
        try:
            return Vector(self.lattice.cartesian_to_fractional(self.r))
        except AttributeError:
            return None

    @rs.setter
    def rs(self, value):
        """Set :math:`x, y, z` components of `Molecule` position vector.

        Parameters
        ----------
        value : array_like
            :math:`x, y, z` coordinates of `Molecule` position vector
            relative to the origin.

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

    def rotate(self, **kwargs):
        """Rotate `Molecule` position vector.

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
        try:
            self.lattice.rotate(**kwargs)
        except AttributeError:
            pass
        super().rotate(**kwargs)

    def translate(self, t, fix_anchor_point=True):
        """Translate `Molecule` position vector by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_point : bool, optional

        """
        if not fix_anchor_point:
            try:
                self.lattice.translate(t)
            except AttributeError:
                pass
        super().translate(t, fix_anchor_point=fix_anchor_point)

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(lattice=self.lattice, xs=self.xs, ys=self.ys,
                               zs=self.zs))
        return super_dict


class LatticeMolecules(Molecules):
    """A `Molecules` sub-class for crystal structure lattice molecules.

    Sub-class of `Molecules` class, and a container class for lists of
    :class:`~sknano.core.molecules.LatticeMolecule` instances.

    Parameters
    ----------
    molecules : {None, sequence, `LatticeMolecules`}, optional
        if not `None`, then a list of `LatticeMolecule` instance objects or an
        existing `LatticeMolecules` instance object.

    """

    @property
    def __molecule_class__(self):
        return LatticeMolecule

    @property
    def rs(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`Molecule.r` position \
            `Vector`\ s"""
        return np.asarray([molecule.rs for molecule in self])

    @property
    def xs(self):
        """:class:`~numpy:numpy.ndarray` of `Molecule`\ s :math:`x` \
            coordinates."""
        return self.rs[:, 0]

    @property
    def ys(self):
        """:class:`~numpy:numpy.ndarray` of `Molecule`\ s :math:`y` \
            coordinates."""
        return self.rs[:, 1]

    @property
    def zs(self):
        """:class:`~numpy:numpy.ndarray` of `Molecule`\ s :math:`z` \
            coordinates."""
        return self.rs[:, 2]

    @property
    def lattice(self):
        try:
            return self[0].lattice
        except IndexError:
            return None

    @lattice.setter
    def lattice(self, value):
        [setattr(molecule, 'lattice', value) for molecule in self]
