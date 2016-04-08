# -*- coding: utf-8 -*-
"""
===============================================================================
Topology Mixin (:mod:`sknano.core.atoms.mixins.topology`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins.topology

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import Iterable, namedtuple
# from operator import attrgetter

# import numpy as np
# np.set_printoptions(edgeitems=20)
# np.set_printoptions(threshold=10000)

from .angles import Angle, Angles, compute_angle
from .bonds import Bond, Bonds, compute_bond
from .dihedrals import Dihedral, Dihedrals, compute_dihedral
from .impropers import Improper, Impropers, compute_improper

__all__ = ['AtomTopologyMixin', 'AtomsTopologyMixin', 'AtomsTopologyStats']

operand_error_msg = 'Expected an `iterable` object containing {}'
ids_operand_error_msg = operand_error_msg.format('`ints`')
AtomsTopologyStats = namedtuple('AtomsTopologyStats',
                                ('angles', 'bonds', 'dihedrals', 'impropers'))


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

    @angles.setter
    def angles(self, value):
        if not isinstance(value, Angles):
            raise TypeError('Expected an `Angles` object')
        self._angles = value

    @property
    def bonds(self):
        """Atom :class:`~sknano.core.atoms.mixins.Bonds`"""
        try:
            return self._bonds
        except AttributeError:
            self._update_bonds()
            return self._bonds

    @bonds.setter
    def bonds(self, value):
        if not isinstance(value, Bonds):
            raise TypeError('Expected a `Bonds` object.')
        self._bonds = value

    @property
    def dihedrals(self):
        """Atom :class:`~sknano.core.atoms.mixins.Dihedrals`"""
        try:
            return self._dihedrals
        except AttributeError:
            self._update_dihedrals()
            return self._dihedrals

    @dihedrals.setter
    def dihedrals(self, value):
        if not isinstance(value, Dihedrals):
            raise TypeError('Expected a `Dihedrals` object')
        self._dihedrals = value

    @property
    def impropers(self):
        """Atom :class:`~sknano.core.atoms.mixins.Dihedrals`"""
        try:
            return self._impropers
        except AttributeError:
            self._update_impropers()
            return self._impropers

    @impropers.setter
    def impropers(self, value):
        if not isinstance(value, Impropers):
            raise TypeError('Expected an `Impropers` object')
        self._impropers = value

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
            bonds = Bonds([Bond(self, nn) for nn in self.NN])
            bonds = bonds.unique
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
    def angles_in_degrees(self):
        """:class:`~python:bool` setting for returning angles in degrees."""
        try:
            return self._angles_in_degrees
        except AttributeError:
            self._angles_in_degrees = False
            return False

    @angles_in_degrees.setter
    def angles_in_degrees(self, value):
        if not isinstance(value, bool):
            raise ValueError('Expected a boolean value.')
        self._angles_in_degrees = value
        self._update_topology()

    @property
    def all_angles(self):
        """:class:`~sknano.core.atoms.mixins.Angles`."""
        try:
            return self._all_angles
        except AttributeError:
            self._update_angles()
            return self._all_angles

    @property
    def all_bonds(self):
        """:class:`~sknano.core.atoms.mixins.Bonds`."""
        try:
            return self._all_bonds
        except AttributeError:
            self._update_bonds()
            return self._all_bonds

    @property
    def all_dihedrals(self):
        """:class:`~sknano.core.atoms.mixins.Dihedrals`."""
        try:
            return self._all_dihedrals
        except AttributeError:
            self._update_dihedrals()
            return self._all_dihedrals

    @property
    def all_impropers(self):
        """:class:`~sknano.core.atoms.mixins.Impropers`."""
        try:
            return self._all_impropers
        except AttributeError:
            self._update_impropers()
            return self._all_impropers

    @property
    def angles(self):
        """:class:`~sknano.core.atoms.mixins.Angles`."""
        return self.all_angles.unique

    @property
    def bonds(self):
        """:class:`~sknano.core.atoms.mixins.Bonds`."""
        return self.all_bonds.unique

    @property
    def dihedrals(self):
        """:class:`~sknano.core.atoms.mixins.Dihedrals`."""
        return self.all_dihedrals.unique

    @property
    def impropers(self):
        """:class:`~sknano.core.atoms.mixins.Impropers`."""
        return self.all_impropers.unique

    @property
    def topology_stats(self):
        """:class:`~python:dict` of topology statistics."""
        topostats = {}
        topostats['angles'] = self.angles.statistics
        topostats['bonds'] = self.bonds.statistics
        topostats['dihedrals'] = self.dihedrals.statistics
        # topostats['impropers'] = self.impropers.statistics
        topostats['impropers'] = None
        return AtomsTopologyStats(**topostats)

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
                bonds = bonds.atom_ids

        elif ids_only:
            bonds = bonds.atom_ids

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

    def _update_topology(self):
        self._update_angles()
        self._update_bonds()
        self._update_dihedrals()
        self._update_impropers()

    def _update_angles(self):
        self._all_angles = \
            Angles([angle for atom in self for angle in atom.angles
                    if len(atom.angles) > 0], parent=self)
        if self.angles_in_degrees:
            self._all_angles.degrees = True

    def _update_bonds(self):
        self._all_bonds = \
            Bonds([bond for atom in self for bond in atom.bonds
                   if len(atom.bonds) > 0], parent=self)

    def _update_dihedrals(self):
        self._all_dihedrals = \
            Dihedrals([dihedral for atom in self for dihedral in
                       atom.dihedrals if len(atom.dihedrals) > 0],
                      parent=self)
        if self.angles_in_degrees:
            self._all_dihedrals.degrees = True

    def _update_impropers(self):
        self._all_impropers = \
            Impropers([improper for atom in self for improper in
                       atom.impropers if len(atom.impropers) > 0],
                      parent=self)
        if self.angles_in_degrees:
            self._all_impropers.degrees = True
