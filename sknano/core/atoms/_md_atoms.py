# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes for MD trajectory analysis (:mod:`sknano.core.atoms._md_atoms`)
===============================================================================

An `Atom` class for molecular dynamics trajectory analysis.

.. currentmodule:: sknano.core.atoms._md_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np

import sknano.core.atoms

from ._energy_atoms import EnergyAtom, EnergyAtoms
from ._force_atoms import ForceAtom, ForceAtoms
from ._structure_atoms import StructureAtom, StructureAtoms
from ._bonds import Bond, Bonds

__all__ = ['MDAtom', 'MDAtoms']


class MDAtom(StructureAtom, ForceAtom, EnergyAtom):
    """An `Atom` class for molecular dynamics trajectory analysis.

    """
    def __init__(self, *args, reference_atom=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.reference_atom = reference_atom

    def __dir__(self):
        attrs = super().__dir__()
        attrs.append('reference_atom')
        return attrs

    @property
    def reference_atom(self):
        return self._reference_atom

    @reference_atom.setter
    def reference_atom(self, value):
        self._reference_atom = value
        if value is not None:
            try:
                self.r0 = self.reference_atom.r
            except AttributeError:
                pass

    @property
    def reference_atom_neighbors(self):
        atoms = []
        neighbors = self.neighbors
        try:
            for atom in self.reference_atom.neighbors:
                if atom.id in neighbors.ids:
                    atoms.append(neighbors[
                                 neighbors.ids.tolist().index(atom.id)])
                else:
                    atomdict = atom.todict()
                    [atomdict.update({k: np.inf}) for k in ('x', 'y', 'z')]
                    atoms.append(self.__class__(**atomdict))
        except AttributeError:
            pass
        return sknano.core.atoms.MDAtoms(atoms)

    @property
    def reference_atom_bonds(self):
        """Return atom `Bonds` instance."""
        try:
            return Bonds([Bond(self, nn) for nn in
                          self.reference_atom_neighbors])
        except (AttributeError, TypeError):
            return Bonds()

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(reference_atom=self.reference_atom,
                               neighbors=self.neighbors))
        return super_dict


class MDAtoms(StructureAtoms, ForceAtoms, EnergyAtoms):
    """An `Atoms` sub-class for molecular dynamics trajectory analysis.

    Sub-class of :class:`~sknano.core.atoms.StructureAtoms` class,
    and a container class for lists of :class:`~sknano.core.atoms.MDAtom`
    instances.

    Parameters
    ----------
    atoms : {None, sequence, `MDAtoms`}, optional
        if not `None`, then a list of `MDAtom` instance objects or an
        existing `MDAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return MDAtom
