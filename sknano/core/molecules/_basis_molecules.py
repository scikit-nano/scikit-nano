# -*- coding: utf-8 -*-
"""
===============================================================================
Basis molecule classes (:mod:`sknano.core.molecules._basis_molecules`)
===============================================================================

.. currentmodule:: sknano.core.molecules._basis_molecules

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._lattice_molecules import LatticeMolecule, LatticeMolecules
from ._periodic_molecules import PBCMolecule, PBCMolecules

__all__ = ['BasisMolecule', 'BasisMolecules']


class BasisMolecule(PBCMolecule, LatticeMolecule):
    """Class representation of a crystal structure basis molecule.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.Crystal3DLattice`
    xs, ys, zs : float
    """
    pass


class BasisMolecules(PBCMolecules, LatticeMolecules):
    """A `Molecules` sub-class for crystal structure basis molecules.

    Sub-class of `Molecules` class, and a container class for lists of
    :class:`~sknano.core.molecules.BasisMolecule` instances.

    Parameters
    ----------
    molecules : {None, sequence, `BasisMolecules`}, optional
        if not `None`, then a list of `BasisMolecule` objects or an
        instance of a `BasisMolecules` object.

    """
    @property
    def __molecule_class__(self):
        return BasisMolecule
