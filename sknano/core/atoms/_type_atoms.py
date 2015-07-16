# -*- coding: utf-8 -*-
"""
===============================================================================
Atoms class for :class:`TypeAtom`\ s (:mod:`sknano.core.atoms._type_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._type_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

# import numpy as np

from sknano.core import dedupe
from ._atoms import Atoms
from ._type_atom import TypeAtom

__all__ = ['TypeAtoms']


class TypeAtoms(Atoms):
    """An eXtended `Atoms` class.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.TypeAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `TypeAtoms`}, optional
        if not `None`, then a list of `TypeAtom` instance objects or an
        existing `TypeAtoms` instance object.

    """
    def __init__(self, atoms=None):

        super().__init__(atoms)
        self._types = {}

    @property
    def __atom_class__(self):
        return TypeAtom

    def sort(self, key=attrgetter('type'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def Ntypes(self):
        """Number of :attr:`~TypeAtoms.types`."""
        return len(list(self.types.keys()))

    @property
    def types(self):
        """:class:`python:dict` of :attr:`TypeAtom.type`\ s."""
        self._update_types()
        return self._types

    @property
    def atomtypes(self):
        return self.types

    def _update_types(self):
        [self.add_type(atom) for atom in self]

    def add_type(self, atom):
        """Add atom type to :attr:`~TypeAtoms.types`.

        Parameters
        ----------
        atom : :class:`~sknano.core.atoms.TypeAtom`
            A :class:`~sknano.core.atoms.TypeAtom` instance.

        """
        self._types[atom.type] = {}
        self._types[atom.type]['mass'] = atom.mass

    def add_atomtype(self, atom):
        self.add_type(atom)

    def add_types(self, atoms=None):
        """Add atom type for each atom in atoms to :attr:`TypeAtom.types` \
            dictionary.

        Parameters
        ----------
        atoms : sequence
            a list of `Atom` object instances

        """
        try:
            [self.add_type(atom) for atom in atoms]
        except TypeError:
            print('Expected an iterable sequence of `Atom` objects.')

    def add_atomtypes(self, atoms=None):
        self.add_types(atoms=atoms)

    def assign_unique_types(self, from_attr='element'):
        """Assign unique :attr:`TypeAtom.type`\s to each `TypeAtom` in \
            `TypeAtoms` from an existing unique atom attribute.

        .. versionchanged:: 0.3.11

           Now accepts a keyword argument `from_attr`.

        The assignment of unique :attr:`TypeAtom.type`\s is performed
        by mapping an existing atom attribute (default: element)
        to a unique integer, starting at 1.

        Parameters
        ----------
        from_attr : :class:`~python:str`
            An existing atom attribute used to generate an attribute
            mapping that maps the attribute to a unique atom
            :attr:`~TypeAtom.type`.

        """
        attrlist = [getattr(atom, from_attr) for atom in self]
        attrmap = \
            {attr: i for i, attr in enumerate(dedupe(attrlist), start=1)}
        self.mapatomattr(from_attr, 'type', attrmap)

    def get_types(self, asdict=False):
        """Return list of `TypeAtom` :attr:`TypeAtom.type`\ s.

        Parameters
        ----------
        asdict : :class:`python:bool`, optional

        Returns
        -------
        :class:`python:list` if `asdict` is `False`
        :class:`python:dict` if `asdict` is `True`

        """
        if asdict:
            return self.types
        else:
            return [atom.type for atom in self]

    def get_atomtypes(self, asdict=False):
        return self.get_types(asdict=asdict)
