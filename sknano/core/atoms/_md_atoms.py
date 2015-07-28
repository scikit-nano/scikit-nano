# -*- coding: utf-8 -*-
"""
===============================================================================
Atoms class for MD trajectory analysis (:mod:`sknano.core.atoms._md_atoms`)
===============================================================================

An `Atoms` class for molecular dynamics trajectory analysis.

.. currentmodule:: sknano.core.atoms._md_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import numbers
# import numpy as np

from ._energy_atoms import EnergyAtoms
from ._force_atoms import ForceAtoms
from ._structure_atoms import StructureAtoms
from ._md_atom import MDAtom

__all__ = ['MDAtoms']


class MDAtoms(StructureAtoms, EnergyAtoms, ForceAtoms):
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
