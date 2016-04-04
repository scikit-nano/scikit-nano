# -*- coding: utf-8 -*-
"""
===============================================================================
Class representations of atom bonds (:mod:`sknano.core.atoms.mixins._bonds`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins._bonds

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import namedtuple

import numpy as np

from sknano.core.math import Vector
from ._topology_base import Topology, TopologyCollection, TopologyStats, \
    check_operands as check_operands_

__all__ = ['compute_bond', 'Bond', 'Bonds', 'BondStats']

BondStats = namedtuple('BondStats', TopologyStats._fields)


def compute_bond(*atoms, check_operands=True, degrees=False):
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
        if len(atoms) != 2.

    """
    if check_operands:
        atoms = check_operands_(*atoms, size=2)
    atom1, atom2 = atoms
    bvec = Bond(atom1, atom2).vector
    return bvec.length


class Bond(Topology):
    """Class representation of bond between 2 `Atom` objects.

    Parameters
    ----------
    origin, end : :class:`~sknano.core.atoms.Atom` instances
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.
    id : :class:`~python:int`
    order : {1, 2, 3, 5}

    """
    def __init__(self, *args, order=None, **kwargs):
        super().__init__(*args, size=2, **kwargs)

        self._order = order
        self._visited = False
        self.fmtstr = "{origin!r}, {end!r}, order={order!r}, " + super().fmtstr

    def __str__(self):
        """Return nice string representation of `Bond`."""
        return "Bond({!r}->{!r})".format(self.origin.id, self.end.id)

    # def __dir__(self):
    #     attrs = super().__dir__()
    #     attrs.extend(['vector', 'unit_vector'])
    #     return attrs
    #     return ['origin', 'end', 'parent', 'vector', 'unit_vector', 'length']
    #     return ['atoms', 'atom1', 'atom2', 'vector', 'unit_vector', 'length']

    @property
    def origin(self):
        """:class:`~sknano.core.atoms.Atom` 1 in `Bond`."""
        return self.atoms[0]

    @property
    def end(self):
        """:class:`~sknano.core.atoms.Atom` 2 in `Bond`."""
        return self.atoms[1]

    @property
    def atom1(self):
        """An alias for :attr:`~Bond.origin`."""
        return self.origin

    @property
    def atom2(self):
        """An alias for :attr:`~Bond.end`."""
        return self.end

    @property
    def vector(self):
        """:class:`Bond` :class:`~sknano.core.math.Vector`.

        :class:`Bond` :class:`~sknano.core.math.Vector` points from
        :attr:`Bond.origin` to :attr:`Bond.end`.

        .. note::
           Accounts for periodic boundary conditions if a
           :class:`~sknano.core.crystallography.Crystal3DLattice` is assigned
           to the :attr:`~Topology.atoms`.

        """
        try:
            # lattice = self.origin.lattice
            if (self.parent is not None and any(self.parent.pbc)) or \
                    self.parent is None:
                lattice = self.atoms.lattice
                dr = lattice.fractional_to_cartesian(
                    lattice.fractional_diff(self.end.rs, self.origin.rs))
        except AttributeError:
            dr = self.end.r - self.origin.r
        return Vector(dr, p0=self.origin.r.p)

    @property
    def unit_vector(self):
        """`Bond` :attr:`~sknano.core.math.Vector.unit_vector`."""
        return self.vector.unit_vector

    @property
    def bond(self):
        """An alias for :attr:`TopologyCollection.measure`."""
        return self.measure

    @property
    def length(self):
        """An alias for :attr:`TopologyCollection.measure`."""
        return self.measure

    @property
    def order(self):
        """Bond order."""
        return self._order

    @order.setter
    def order(self, value):
        self._order = value

    def compute_measure(self):
        """`Bond` :attr:`~sknano.core.math.Vector.length`."""
        return self.vector.length

    def partner(self, atom):
        """Return :class:`~sknano.core.atoms.Atom` bonded to `atom`."""
        if atom not in self.atoms:
            return None
        if atom is self.origin:
            return self.end
        else:
            return self.origin

    # def ring_groups(self):
    #     return self.rings

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        super_dict = super().todict()
        super_dict.update(dict(origin=self.origin, end=self.end,
                               order=self.order))
        return super_dict


class Bonds(TopologyCollection):
    """Base class for collection of atom `Bond`\ s.

    Parameters
    ----------
    topolist : {None, sequence, `Bonds`}, optional
        if not `None`, then a list of `Bond` objects
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.

    """
    @property
    def __item_class__(self):
        return Bond

    @property
    def Nbonds(self):
        """Number of `Bond`\ s in `Bonds`."""
        return len(self)

    @property
    def bonds(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`~Bond.bond`\ s."""
        return self.measures

    @property
    def lengths(self):
        """An alias for :attr:`Bonds.bonds`."""
        return self.measures

    @property
    def mean_bond(self):
        """Mean bond length."""
        return self.mean_measure

    @property
    def mean_length(self):
        """An alias for :attr:`Bonds.mean_bond`."""
        return self.mean_bond

    @property
    def vectors(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`~Bond.vector`\ s."""
        return np.asarray([bond.vector for bond in self])

    @property
    def unit_vectors(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`~Bond.unit_vector`\ s."""
        return np.asarray([bond.unit_vector for bond in self])

    @property
    def statistics(self):
        """Bond stats."""
        return BondStats(**super().statistics._asdict())
