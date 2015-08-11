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

__all__ = ['MDAtom', 'MDAtoms']


class MDAtom(StructureAtom, ForceAtom, EnergyAtom):
    """An `Atom` class for molecular dynamics trajectory analysis.

    """
    def __init__(self, *args, reference_atom=None, t0_atom=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.reference_atom = reference_atom
        self.t0_atom = t0_atom
        try:
            self.r0 = t0_atom.r
        except AttributeError:
            pass

    def __dir__(self):
        attrs = super().__dir__()
        # attrs.append('reference_atom')
        attrs.extend(['reference_atom', 't0_atom'])
        return attrs

    @property
    def NN(self):
        return super().NN

    @NN.setter
    def NN(self, value):
        """Set nearest-neighbor `Atoms`."""
        try:
            atoms = []
            for atom in self.reference_atom.NN:
                if atom.id in value.ids:
                    atoms.append(value[value.ids.tolist().index(atom.id)])
                else:
                    atomdict = atom.todict()
                    [atomdict.update({k: np.inf}) for k in ('x', 'y', 'z')]
                    print(atomdict)
                    atoms.append(self.__class__(**atomdict))
            value = sknano.core.atoms.MDAtoms(atoms)
        except AttributeError:
            pass
        super(MDAtom, MDAtom).NN.__set__(self, value)

    @NN.deleter
    def NN(self):
        super(MDAtom, MDAtom).NN.__delete__(self)

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(reference_atom=self.reference_atom,
                               t0_atom=self.t0_atom,
                               NN=self.NN))
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
