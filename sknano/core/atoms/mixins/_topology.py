# -*- coding: utf-8 -*-
"""
===============================================================================
Topology Mixin (:mod:`sknano.core.atoms.mixins._topology`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins._topology

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import Iterable
# from operator import attrgetter

# import numpy as np
# np.set_printoptions(edgeitems=20)
# np.set_printoptions(threshold=10000)

from ._angles import Angle, Angles, compute_angle
from ._bonds import Bond, Bonds, compute_bond
from ._dihedrals import Dihedral, Dihedrals, compute_dihedral
from ._impropers import Improper, Impropers, compute_improper

__all__ = ['AtomTopologyMixin', 'AtomsTopologyMixin']

operand_error_msg = 'Expected an `iterable` object containing {}'
ids_operand_error_msg = operand_error_msg.format('`ints`')


class AtomTopologyMixin:
    """Mixin :class:`~sknano.core.atoms.Atom` topology class."""
    @property
    def angles(self):
        """Atom bond :class:`~sknano.core.atoms.mixins.Angles`."""
        try:
            return self._angles
        except AttributeError:
            self._update_angles()
            return self._angles

    @property
    def bonds(self):
        """Atom :class:`~sknano.core.atoms.mixins.Bonds`"""
        try:
            return self._bonds
        except AttributeError:
            self._update_bonds()
            return self._bonds

    @property
    def dihedrals(self):
        """Atom :class:`~sknano.core.atoms.mixins.Dihedrals`"""
        try:
            return self._dihedrals
        except AttributeError:
            self._update_dihedrals()
            return self._dihedrals

    @property
    def impropers(self):
        """Atom :class:`~sknano.core.atoms.mixins.Dihedrals`"""
        try:
            return self._impropers
        except AttributeError:
            self._update_impropers()
            return self._impropers

    def _update_angles(self):
        try:
            angles = Angles()
            neighbors = [bond.partner(self) for bond in self.bonds]
            for lneighbor in neighbors:
                for rneighbor in neighbors:
                    if lneighbor != rneighbor:
                        angles.append(Angle(lneighbor, self, rneighbor))
            angles = angles.unique
        except (AttributeError, TypeError):
            angles = Angles()
        self._angles = angles

    def _update_bonds(self):
        try:
            bonds = Bonds([Bond(self, nn) for nn in self.NN]).unique
        except (AttributeError, TypeError):
            bonds = Bonds()
        self._bonds = bonds

    def _update_dihedrals(self):
        try:
            dihedrals = Dihedrals()
            bonds = self.bonds
            for bond in bonds:
                a1, a2 = bond.atoms
                a1neighbors = [bond.partner(a1) for bond in a1.bonds]
                a2neighbors = [bond.partner(a2) for bond in a2.bonds]
                for a1n in a1neighbors:
                    if a1n not in (a1, a2):
                        for a2n in a2neighbors:
                            if a2n not in (a1, a2) and a1n != a2n:
                                dihedrals.append(Dihedral(a1n, a1, a2, a2n))
            dihedrals = dihedrals.unique
        except (AttributeError, TypeError):
            dihedrals = Dihedrals()
        self._dihedrals = dihedrals

    def _update_impropers(self):
        try:
            impropers = Impropers()
            for dihedral in self.dihedrals:
                impropers.append(Improper(*dihedral.atoms))
            impropers = impropers.unique
        except (AttributeError, TypeError):
            impropers = Impropers()
        self._impropers = impropers


class AtomsTopologyMixin:
    """Mixin :class:`~sknano.core.atoms.Atoms` topology class."""
    @property
    def angles(self):
        """:class:`~sknano.core.atoms.mixins.Angles`."""
        return Angles([angle for atom in self for angle in atom.angles
                       if len(atom.angles) > 0], parent=self)

    @property
    def bonds(self):
        """:class:`~sknano.core.atoms.mixins.Bonds`."""
        return Bonds([bond for atom in self for bond in atom.bonds
                      if len(atom.bonds) > 0], parent=self)

    @property
    def dihedrals(self):
        """:class:`~sknano.core.atoms.mixins.Dihedrals`."""
        return Dihedrals([dihedral for atom in self for dihedral in
                          atom.dihedrals if len(atom.dihedrals) > 0],
                         parent=self)

    @property
    def impropers(self):
        """:class:`~sknano.core.atoms.mixins.Impropers`."""
        return Impropers([improper for atom in self for improper in
                          atom.impropers if len(atom.impropers) > 0],
                         parent=self)

    @property
    def topology_stats(self):
        """:class:`~python:dict` of topology statistics."""
        stats = {}
        stats['Angles'] = self.angles.statistics
        stats['Bonds'] = self.bonds.statistics
        stats['Dihedrals'] = self.dihedrals.statistics
        stats['Impropers'] = self.impropers.statistics
        return stats

    def get_angle(self, *triplet, check_operands=True, degrees=False):
        """Compute bond angles.

        Parameters
        ----------
        triplet : :class:`~python:list` of :class:`~python:int`\ s.
            :class:`~python:list` of :attr:`~sknano.core.atoms.IDAtom.id`\ s
        check_operands : :class:`~python:bool`, optional
        degrees : :class:`~python:bool`, optional

        Returns
        -------
        :class:`~python:float`

        Raises
        ------
        :class:`~python:TypeError`
            if `triplet` is not a list of
            :attr:`~sknano.core.atoms.IDAtom.id`\ s
        :class:`~python:ValueError`
            if len(triplet) != 3.

        """
        if check_operands:
            if not isinstance(triplet, Iterable):
                raise TypeError(ids_operand_error_msg)
            if len(triplet) == 1:
                if isinstance(triplet[0], Iterable):
                    return self.get_angle(*triplet[0], degrees=degrees)
                else:
                    triplet = triplet[0]

            if not isinstance(triplet, Iterable) or not \
                    all([isinstance(id, int) for id in triplet]):
                raise TypeError(ids_operand_error_msg)

            if len(triplet) != 3:
                raise ValueError('Expected 3 atom ids for bond angle '
                                 'calculation')

        return compute_angle(*self.get_atoms(ids=triplet), degrees=degrees,
                             check_operands=False)

    def get_bond(self, *pair, check_operands=True, degrees=False):
        """Compute bond lengths.

        Parameters
        ----------
        pair : :class:`~python:list` of :class:`~python:int`\ s.
            :class:`~python:list` of :attr:`~sknano.core.atoms.IDAtom.id`\ s
        check_operands : :class:`~python:bool`, optional
        degrees : :class:`~python:bool`, optional

        Returns
        -------
        :class:`~python:float`

        Raises
        ------
        :class:`~python:TypeError`
            if `pair` is not a list of
            :attr:`~sknano.core.atoms.IDAtom.id`\ s
        :class:`~python:ValueError`
            if len(pair) != 2.

        """
        if check_operands:
            if not isinstance(pair, Iterable):
                raise TypeError(ids_operand_error_msg)
            if len(pair) == 1:
                if isinstance(pair[0], Iterable):
                    return self.get_bond(*pair[0], degrees=degrees)
                else:
                    pair = pair[0]

            if not isinstance(pair, Iterable) or not \
                    all([isinstance(id, int) for id in pair]):
                raise TypeError(ids_operand_error_msg)

            if len(pair) != 2:
                raise ValueError('Expected 2 atom ids for bond length '
                                 'calculation')

        return compute_bond(*self.get_atoms(ids=pair), degrees=degrees,
                            check_operands=False)

    def get_dihedral(self, *ids, check_operands=True, degrees=False):
        """Compute dihedral angles."""
        if check_operands:
            if not isinstance(ids, Iterable):
                raise TypeError(ids_operand_error_msg)
            if len(ids) == 1:
                if isinstance(ids[0], Iterable):
                    return self.get_angle(*ids[0], degrees=degrees)
                else:
                    ids = ids[0]

            if not isinstance(ids, Iterable) or not \
                    all([isinstance(id, int) for id in ids]):
                raise TypeError(ids_operand_error_msg)

            if len(ids) != 4:
                raise ValueError('Expected 4 atom ids for dihedral '
                                 'calculation')

        return compute_dihedral(*self.get_atoms(ids=ids), degrees=degrees,
                                check_operands=False)

    def get_improper(self, *ids, check_operands=True, degrees=False):
        """Compute improper angles."""
        if check_operands:
            if not isinstance(ids, Iterable):
                raise TypeError(ids_operand_error_msg)
            if len(ids) == 1:
                if isinstance(ids[0], Iterable):
                    return self.get_angle(*ids[0], degrees=degrees)
                else:
                    ids = ids[0]

            if not isinstance(ids, Iterable) or not \
                    all([isinstance(id, int) for id in ids]):
                raise TypeError(ids_operand_error_msg)

            if len(ids) != 4:
                raise ValueError('Expected 4 atoms for improper calculation')

        return compute_improper(*self.get_atoms(ids=ids), degrees=degrees,
                                check_operands=False)

    def get_angles(self, neighbors_only=False, unique=False, ids_only=False):
        """Return list of bond angles."""
        angles = self.angles

        if isinstance(angles, Angles):
            angles = angles.data

        return angles

    def get_bonds(self, neighbors_only=False, unique=False, ids_only=False):
        """Return list of bonds."""
        bonds = self.bonds
        if neighbors_only and not unique:
            # bonds = [atom.neighbors.data for atom in self]
            bonds = [[bond.partner(atom) for bond in atom.bonds]
                     for atom in self]
            if ids_only:
                bonds = [[neighbor.id for neighbor in neighborlist]
                         for neighborlist in bonds]

        elif unique:
            bonds = bonds.unique
            if neighbors_only and ids_only:
                bonds = [bond[-1].id for bond in bonds]
            elif neighbors_only:
                bonds = [bond[-1] for bond in bonds]
            elif ids_only:
                bonds = bonds.ids

        elif ids_only:
            bonds = bonds.ids

        if isinstance(bonds, Bonds):
            bonds = bonds.data
        return bonds

    def get_dihedrals(self, neighbors_only=False, unique=False,
                      ids_only=False):
        """Return list of :class:`~sknano.core.atoms.mixins.Dihedral`\ s"""
        dihedrals = self.dihedrals

        if isinstance(dihedrals, Dihedrals):
            dihedrals = dihedrals.data

        return dihedrals

    def get_impropers(self, neighbors_only=False, unique=False,
                      ids_only=False):
        """Return list of :class:`~sknano.core.atoms.mixins.Impropers`\ s"""
        impropers = self.impropers

        if isinstance(impropers, Impropers):
            impropers = impropers.data

        return impropers

    def update_attrs(self, topology=True, **kwargs):
        """Update :class:`AtomsTopologyMixin` class attributes."""
        super().update_attrs(**kwargs)
        if topology:
            [atom._update_bonds() for atom in self]
            [atom._update_angles() for atom in self]
            [atom._update_dihedrals() for atom in self]
            [atom._update_impropers() for atom in self]
