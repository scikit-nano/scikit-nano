# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class for MD trajectory analysis (:mod:`sknano.core.atoms._md_atom`)
===============================================================================

An `Atom` class for molecular dynamics trajectory analysis.

.. currentmodule:: sknano.core.atoms._md_atom

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np

import sknano.core.atoms

from ._structure_atoms import StructureAtom as Atom

__all__ = ['MDAtom']


class MDAtom(Atom):
    """An `Atom` class for molecular dynamics trajectory analysis.

    """
    def __init__(self, reference_atom=None, **kwargs):
        super(MDAtom, self).__init__(**kwargs)
        self.reference_atom = reference_atom

    @property
    def NN(self):
        """Nearest-neighbor `Atoms`."""
        if self.reference_atom is not None:
            NN = self._NN
            try:
                if not np.allclose(NN.atom_ids,
                                   self.reference_atom.NN.atom_ids):
                    self.NN = sknano.core.atoms.MDAtoms(
                        atoms=[NN[NN.atom_ids.tolist().index(atom_id)] for
                               atom_id in self.reference_atom.NN.atom_ids])
            except ValueError:
                pass

        return self._NN

    @NN.setter
    def NN(self, value):
        """Set nearest-neighbor `Atoms`."""
        if not isinstance(value, sknano.core.atoms.Atoms):
            raise TypeError('Expected an `Atoms` object.')
        self._NN = value
