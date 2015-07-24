# -*- coding: utf-8 -*-
"""
===============================================================================
`EnergyAtom` container class (:mod:`sknano.core.atoms._energy_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._energy_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from ._atoms import Atoms
from ._energy_atom import EnergyAtom

__all__ = ['EnergyAtoms']


class EnergyAtoms(Atoms):
    """An `Atoms` sub-class for `EnergyAtom`\ s.

    A container class for class:`~sknano.core.atoms.EnergyAtom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `EnergyAtoms`}, optional
        if not `None`, then a list of `EnergyAtom` instance objects or an
        existing `EnergyAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return EnergyAtom

    def sort(self, key=attrgetter('etotal'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def kinetic_energies(self):
        """:class:`~numpy:numpy.ndarray` of `EnergyAtom.ke`."""
        return np.asarray([atom.ke for atom in self])

    @property
    def potential_energies(self):
        """:class:`~numpy:numpy.ndarray` of `EnergyAtom.pe`."""
        return np.asarray([atom.pe for atom in self])

    @property
    def total_energies(self):
        """:class:`~numpy:numpy.ndarray` of `EnergyAtom.etotal`."""
        return np.asarray([atom.etotal for atom in self])
