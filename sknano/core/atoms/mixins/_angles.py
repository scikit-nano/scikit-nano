# -*- coding: utf-8 -*-
"""
===============================================================================
Class representations of bond angles (:mod:`sknano.core.atoms.mixins._angles`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins._angles

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.core.math import vector as vec
from ._bonds import Bond, Bonds
from ._topology_base import AtomTopology, AtomsTopology, \
    check_operands as check_operands_

__all__ = ['compute_angle', 'Angle', 'Angles']


def compute_angle(*atoms, check_operands=True, degrees=False):
    """Compute bond angles.

    Parameters
    ----------
    *atoms : {:class:`~python:list`, :class:`~sknano.core.atoms.Atoms`}
        :class:`~python:list` of :class:`~sknano.core.atoms.Atom`\ s
        or an :class:`~sknano.core.atoms.Atoms` object.
    check_operands : :class:`~python:bool`, optional
    degrees : :class:`~python:bool`, optional

    Returns
    -------
    :class:`~python:float`

    Raises
    ------
    :class:`~python:TypeError`
        if `atoms` is not a list of :class:`~sknano.core.atoms.Atom` objects
        or an :class:`~sknano.core.atoms.Atoms` object.
    :class:`~python:ValueError`
        if len(atoms) != 3.

    """
    if check_operands:
        atoms = check_operands_(*atoms, size=3)

    atom1, atom2, atom3 = atoms
    b21 = Bond(atom2, atom1).vector
    b23 = Bond(atom2, atom3).vector
    return vec.angle(b21, b23, degrees=degrees)


class Angle(AtomTopology):
    """Class representation of bond angle between 3 `Atom` objects.

    Parameters
    ----------
    *atoms : {:class:`~python:list`, :class:`~sknano.core.atoms.Atoms`}
        :class:`~python:list` of :class:`~sknano.core.atoms.Atom`\ s
        or an :class:`~sknano.core.atoms.Atoms` object.
    check_operands : :class:`~python:bool`, optional
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.
    id : :class:`~python:int`

    Raises
    ------
    :class:`~python:TypeError`
        if `atoms` is not a list of :class:`~sknano.core.atoms.Atom` objects
        or an :class:`~sknano.core.atoms.Atoms` object.
    :class:`~python:ValueError`
        if len(atoms) != 3.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if self.check_operands and len(self.atoms) != 3:
                raise ValueError('Expected 3 atoms')

        self.fmtstr = "{lneighbor!r}, {center!r}, {rneighbor!r}, " + \
            super().fmtstr

    def __str__(self):
        """Return nice string representation of `Angle`."""
        return "Angle({!r}<-{!r}->{!r})".format(self.lneighbor.id,
                                                self.center.id,
                                                self.rneighbor.id)

    @property
    def lneighbor(self):
        """:class:`~sknano.core.atoms.Atom` 1 in `Angle`."""
        return self.atoms[0]

    @property
    def center(self):
        """:class:`~sknano.core.atoms.Atom` 2 in `Angle`."""
        return self.atoms[1]

    @property
    def rneighbor(self):
        """:class:`~sknano.core.atoms.Atom` 3 in `Angle`."""
        return self.atoms[2]

    @property
    def atom1(self):
        """An alias for :attr:`~Angle.lneighbor`."""
        return self.lneighbor

    @property
    def atom2(self):
        """An alias for :attr:`~Angle.center`."""
        return self.center

    @property
    def atom3(self):
        """An alias for :attr:`~Angle.rneighbor`."""
        return self.rneighbor

    @property
    def lbond(self):
        """:class:`Bond` from :attr:`~Angle.center` -> \
            :attr:`~Angle.lneighbor`."""
        return Bond(self.center, self.lneighbor)

    @property
    def rbond(self):
        """:class:`Bond` from :attr:`~Angle.center` -> \
            :attr:`~Angle.rneighbor`."""
        return Bond(self.center, self.rneighbor)

    @property
    def angle(self):
        """An alias for :attr:`AtomTopology.measure`."""
        return self.measure

    # @property
    # def vector(self):
    #     """:class:`Angle` :class:`~sknano.core.math.Vector`.

    #     :class:`Angle` :class:`~sknano.core.math.Vector` points from
    #     :attr:`Angle.origin` to :attr:`Angle.end`.

    #     .. note::
    #        Accounts for periodic boundary conditions if a
    #        :class:`~sknano.core.crystallography.Crystal3DLattice` is assigned
    #        to the :attr:`~Angle.atoms`.

    #     """
    #     try:
    #         # lattice = self.origin.lattice
    #         if any(self.pbc):
    #             lattice = self.atoms.lattice
    #             dr = lattice.fractional_to_cartesian(
    #                 lattice.fractional_diff(self.end.rs, self.origin.rs))
    #     except AttributeError:
    #         dr = self.end.r - self.origin.r
    #     return Vector(dr, p0=self.origin.r.p)

    def compute_measure(self):
        """Compute the bond angle, which is the measure of an :class:`Angle`.

        Returns
        -------
        :class:`~python:float`

        """
        return compute_angle(*self.atoms)

    @property
    def bond_pairs(self):
        """`cyclic_pairs` of `Bond`\ s."""
        return (self.lbond, self.rbond)

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        super_dict = super().todict()
        super_dict.update(dict(lneighbor=self.lneighbor, center=self.center,
                               rneighbor=self.rneighbor))
        return super_dict


class Angles(AtomsTopology):
    """`AtomsTopology` sub-class for collection of atom `Angle`\ s.

    Parameters
    ----------
    topolist : {None, sequence, `Angles`}, optional
        if not `None`, then a list of `Angle` objects
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.

    """
    @property
    def __item_class__(self):
        return Angle

    @property
    def Nangles(self):
        """Number of `Angle`\ s in `Angles`."""
        return len(self)

    @property
    def angles(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`~Angle.angle`\ s."""
        return self.measures

    @property
    def mean_angle(self):
        """Mean bond angle."""
        return self.mean_measure

    @property
    def bond_pairs(self):
        """:class:`~python:list` of :attr:`Angle.bond_pairs`"""
        return [Bonds(list(angle.bond_pairs)) for angle in self]
