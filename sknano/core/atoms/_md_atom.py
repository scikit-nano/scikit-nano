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
        super().__init__(**kwargs)
        self.reference_atom = reference_atom

    @property
    def NN(self):
        """Nearest-neighbor `Atoms`."""
        return self._NN

    @NN.setter
    def NN(self, value):
        """Set nearest-neighbor `Atoms`."""
        if not isinstance(value, sknano.core.atoms.Atoms):
            raise TypeError('Expected an `Atoms` object.')
        self._NN = value

        if self.reference_atom is not None:
            self._NN = sknano.core.atoms.MDAtoms(
                atoms=[self.NN[self.NN.ids.tolist().index(atom.id)]
                       if atom.id in self.NN.ids else
                       self.__class__(reference_atom=atom.reference_atom,
                                      element=atom.element, id=atom.id,
                                      mol=atom.mol, type=atom.type,
                                      mass=atom.mass, q=atom.q,
                                      x=np.inf, y=np.inf, z=np.inf)
                       for atom in self.reference_atom.NN])
