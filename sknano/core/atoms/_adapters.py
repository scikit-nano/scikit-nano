# -*- coding: utf-8 -*-
"""
===============================================================================
Atom adapter mixins for 3rd-party compat (:mod:`sknano.core.atoms._adapters`)
===============================================================================

.. currentmodule:: sknano.core.atoms._adapters

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

# import numbers

import numpy as np

# from ._atoms import Atom, Atoms

__all__ = ['AtomAdapterMixin', 'AtomsAdapterMixin',
           'VMDAtomAdapterMixin', 'VMDAtomsAdapterMixin']


class VMDAtomAdapterMixin:
    """An `Atom` mixin adapter class with atom attributes consistent
    with VMD atom attributes.

    """

    @property
    def vmd_atom_id(self):
        """VMD atom index attribute."""
        return self.id - 1

    @vmd_atom_id.setter
    def vmd_atom_id(self, value):
        """Set VMD atom index attribute.

        Parameters
        ----------
        value : int
            atom index

        """
        self.id = value + 1


class VMDAtomsAdapterMixin:
    """An `Atoms` mixin adapter class with atom attributes consistent
    with VMD atom attributes."""

    @property
    def vmd_atom_ids(self):
        """Array of VMD atom index attribute."""
        return np.asarray(self.ids) - 1

    def filter_vmd_atom_ids(self, vmd_atom_ids, invert=False):
        """Filter `Atoms` by :attr:`VMDAtomsAdapterMixin.ids` in
        `vmd_atom_ids`.

        Parameters
        ----------
        vmd_atom_ids : array_like
        invert : bool, optional

        """
        mask = \
            np.in1d(self.vmd_atom_ids, vmd_atom_ids, invert=invert).nonzero()
        self.data = np.asarray(self)[mask].tolist()

    def filtered_vmd_atom_ids(self, vmd_atom_ids, invert=False):
        """Return new `Atoms` filtered by :attr:`VMDAtomsAdapterMixin.ids` in
        `vmd_atom_ids`.

        Parameters
        ----------
        vmd_atom_ids : array_like
        invert : bool, optional

        Returns
        -------
        filtered_atoms : `Atoms`
            An instance of `Atoms` (sub)class.

        """
        mask = \
            np.in1d(self.vmd_atom_ids, vmd_atom_ids, invert=invert).nonzero()
        return self.__class__(atoms=np.asarray(self)[mask].tolist(),
                              **self.kwargs)

    def get_atom_from_vmd_atom_id(self, vmd_atom_id):
        """Get `Atom` by VMD atom index.

        Parameters
        ----------
        vmd_atom_id : :class:`~python:int`

        Returns
        -------
        :class:`~sknano.core.atoms.Atom` or `None`

        """
        try:
            return self[np.where(self.vmd_atom_ids == vmd_atom_id)[0]]
        except TypeError:
            print('No atom with id = {}'.format(vmd_atom_id))
            return None


class AtomAdapterMixin(VMDAtomAdapterMixin):
    pass


class AtomsAdapterMixin(VMDAtomsAdapterMixin):
    pass
