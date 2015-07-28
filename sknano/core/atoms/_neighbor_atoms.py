# -*- coding: utf-8 -*-
"""
===============================================================================
NeighborAtom container class (:mod:`sknano.core.atoms._neighbor_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._neighbor_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import numpy as np

from ._cn_atoms import CNAtoms
from ._id_atoms import IDAtoms
from ._xyz_atoms import XYZAtoms
# from ._bonds import Bonds
# from ._neighbor_atom import NeighborAtom
from ._structure_atom import StructureAtom

__all__ = ['NeighborAtoms']


class NeighborAtoms(CNAtoms, IDAtoms, XYZAtoms):
    """An `Atoms` sub-class for nearest-neighbor structure analysis.

    Parameters
    ----------
    atoms : {None, sequence, `NeighborAtoms`}, optional
        if not `None`, then a list of `NeighborAtom` instance objects or an
        existing `NeighborAtoms` instance object.
    kNN : :class:`~python:int`
        Number of nearest neighbors to return when querying the kd-tree.
    NNrc : :class:`~python:float`
        Nearest neighbor radius cutoff.

    """
    @property
    def __atom_class__(self):
        # return NeighborAtom
        return StructureAtom
