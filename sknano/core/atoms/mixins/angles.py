# -*- coding: utf-8 -*-
"""
===============================================================================
Class representations of bond angles (:mod:`sknano.core.atoms.mixins.angles`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins.angles

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import namedtuple

from sknano.core.math import vector as vec
from .bonds import Bond, Bonds
from .topology_base import AngularTopology, AngularTopologyCollection, \
    TopologyStats, check_operands as check_operands_

__all__ = ['compute_angle', 'Angle', 'Angles', 'AngleStats']

AngleStats = namedtuple('AngleStats', TopologyStats._fields)


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


class Angle(AngularTopology):
    """Class representation of bond angle between 3 `Atom` objects.

    Parameters
    ----------
    *atoms : {:class:`~python:list`, :class:`~sknano.core.atoms.Atoms`}
        :class:`~python:list` of :class:`~sknano.core.atoms.Atom`\ s
        or an :class:`~sknano.core.atoms.Atoms` object.
    size : :class:`~python:int`
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.
    id : :class:`~python:int`
    check_operands : :class:`~python:bool`, optional
    degrees : :class:`~python:bool`, optional

    Raises
    ------
    :class:`~python:TypeError`
        if `atoms` is not a list of :class:`~sknano.core.atoms.Atom` objects
        or an :class:`~sknano.core.atoms.Atoms` object.
    :class:`~python:ValueError`
        if len(atoms) != 3.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, size=3, **kwargs)

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

    def compute_measure(self):
        """Compute the bond angle, which is the measure of an :class:`Angle`.

        Returns
        -------
        :class:`~python:float`

        """
        return compute_angle(*self.atoms, degrees=self.degrees)

    @property
    def bond_pairs(self):
        """:class:`~python:tuple` of `Angle` \
            :attr:`Angle.lbond`, :attr:`Angle.rbond`."""
        return (self.lbond, self.rbond)

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        super_dict = super().todict()
        super_dict.update(dict(lneighbor=self.lneighbor, center=self.center,
                               rneighbor=self.rneighbor))
        return super_dict


class Angles(AngularTopologyCollection):
    """`TopologyCollection` sub-class for collection of atom `Angle`\ s.

    Parameters
    ----------
    topolist : {None, sequence, `Angles`}, optional
        if not `None`, then a list of `Angle` objects
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.
    degrees : :class:`~python:bool`, optional

    """
    @property
    def __item_class__(self):
        return Angle

    @property
    def Nangles(self):
        """Number of `Angle`\ s in `Angles`."""
        return len(self)

    @property
    def bond_pairs(self):
        """:class:`~python:list` of :attr:`Angle.bond_pairs`"""
        return [Bonds(list(angle.bond_pairs)) for angle in self]

    @property
    def statistics(self):
        """Angle stats."""
        return AngleStats(**super().statistics._asdict())
