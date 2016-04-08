# -*- coding: utf-8 -*-
"""
===============================================================================
Mixins for 3rd-party compat. (:mod:`sknano.core.atoms.mixins.adapters`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins.adapters

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

# import numbers

import numpy as np

# from .atoms import Atom, Atoms

__all__ = ['AtomAdapterMixin', 'AtomsAdapterMixin',
           'VMDAtomAdapterMixin', 'VMDAtomsAdapterMixin']


vmd_keyword_map = {}
vmd_keyword_map['index'] = 'vmd_indices'
vmd_keyword_map['serial'] = 'ids'
vmd_keyword_map['type'] = 'types'
vmd_keyword_map['element'] = 'elements'


class VMDAtomAdapterMixin:
    """An `Atom` mixin adapter class with atom attributes consistent
    with VMD atom attributes.

    """

    @property
    def vmd_index(self):
        """VMD atom index attribute."""
        return self.id - 1


class VMDAtomsAdapterMixin:
    """An `Atoms` mixin adapter class with atom attributes consistent
    with VMD atom attributes.


    Notes
    -----
    The VMD indices are only consistent with VMD when the number of atoms
    equals the maximum :class:`~sknano.core.atoms.IDAtom.id`.

    """

    @property
    def vmd_indices(self):
        """Array of VMD atom index attribute."""
        return np.asarray(self.ids) - 1

    def filter_vmd_indices(self, vmd_indices, invert=False):
        """Filter `Atoms` by :attr:`VMDAtomsAdapterMixin.vmd_indices`.

        Parameters
        ----------
        vmd_indices : array_like
        invert : bool, optional

        """
        mask = \
            np.in1d(self.vmd_indices, vmd_indices, invert=invert).nonzero()
        self.data = np.asarray(self)[mask].tolist()

    def filtered_vmd_indices(self, vmd_indices, invert=False):
        """Return new `Atoms` filtered by `vmd_indices`.

        Parameters
        ----------
        vmd_indices : array_like
        invert : bool, optional

        Returns
        -------
        filtered_atoms : `Atoms`
            An instance of `Atoms` (sub)class.

        """
        mask = \
            np.in1d(self.vmd_indices, vmd_indices, invert=invert).nonzero()
        return self.__class__(atoms=np.asarray(self)[mask].tolist(),
                              **self.kwargs)

    def get_atom_from_vmd_index(self, vmd_index):
        """Get `Atom` by VMD atom index.

        Parameters
        ----------
        vmd_index : :class:`~python:int`

        Returns
        -------
        :class:`~sknano.core.atoms.Atom` or `None`

        """
        try:
            return self[np.where(self.vmd_indices == vmd_index)[0]]
        except TypeError:
            print('No atom with id = {}'.format(vmd_index))
            return None

    def get_vmd_selection_string(self, keyword):
        """Get a VMD selection string for the VMD keyword."""
        attr = vmd_keyword_map.get(keyword, keyword)
        if hasattr(self, attr):
            return ' '.join((keyword, ' '.join(map(str, getattr(self, attr)))))
        return None


class AtomAdapterMixin(VMDAtomAdapterMixin):
    """Mixin `Atom` class for 3rd party package compatibility."""
    pass


class AtomsAdapterMixin(VMDAtomsAdapterMixin):
    """Mixin `Atoms` class for 3rd party package compatibility."""
    pass
