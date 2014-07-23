# -*- coding: utf-8 -*-
"""
==============================================================================
Class for converting Atoms (:mod:`sknano.structure_io.atoms._atoms_converter`)
==============================================================================

.. currentmodule:: sknano.structure_io.atoms._atoms_converter

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from inspect import getargspec

from ._lammps_atom import LAMMPSAtom
from ._lammps_atoms import LAMMPSAtoms
from ._xyz_atom import XYZAtom
from ._xyz_atoms import XYZAtoms

__all__ = ['AtomsConverter']


class AtomsConverter(object):
    """Class for converting `Atom` objects from one type to another.

    Parameters
    ----------
    atoms : sequence
        Atoms object.
    to : str
        Name of output atoms type

    """

    def __init__(self, atoms=None, to=None):
        self._atoms = None
        Atom = None

        if to.lower() == 'lammps':
            self._atoms = LAMMPSAtoms()
            Atom = LAMMPSAtom
        elif to.lower() == 'xyz':
            self._atoms = XYZAtoms()
            Atom = XYZAtom

        atom_argspec = getargspec(Atom.__init__)
        atom_argdict = dict(zip(atom_argspec.args[1:], atom_argspec.defaults))

        for atom in atoms:
            atom_args = atom_argdict.copy()
            for k, v in atom.atomdict.iteritems():
                if k in atom_args:
                    atom_args[k] = v
            self._atoms.append(Atom(**atom_args))
            del atom_args

    @property
    def atoms(self):
        """Return converted atoms."""
        return self._atoms
