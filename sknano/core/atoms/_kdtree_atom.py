# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class for KDTree analysis (:mod:`sknano.core.atoms._kdtree_atom`)
===============================================================================

.. currentmodule:: sknano.core.atoms._kdtree_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np

import sknano.core.atoms
from ._cn_atom import CNAtom
from ._id_atom import IDAtom
from ._xyz_atom import XYZAtom
from ._bond import Bond
from ._bonds import Bonds

__all__ = ['KDTAtom']


class KDTAtom(CNAtom, XYZAtom, IDAtom):
    """An `Atom` class for KDTree analysis."""
    def __init__(self, *args, NN=None, **kwargs):
        super().__init__(*args, **kwargs)

        self._neighbors = None
        if NN is not None:
            self.NN = NN
        # self.fmtstr = super().fmtstr + ", NN={NN!r}"

    def __dir__(self):
        attrs = super().__dir__()
        attrs.append('NN')
        # attrs.extend(['NN', 'bonds'])
        return attrs

    @CNAtom.CN.getter
    def CN(self):
        """`KDTAtom` coordination number."""
        try:
            return self.NN.Natoms
        except AttributeError:
            return super().CN

    @property
    def NN(self):
        """Nearest-neighbor `Atoms`."""
        try:
            return self._NN
        except AttributeError:
            return None

    @NN.setter
    def NN(self, value):
        """Set nearest-neighbor `Atoms`."""
        if not isinstance(value, sknano.core.atoms.Atoms):
            raise TypeError('Expected an `Atoms` object.')
        self._NN = value

    @NN.deleter
    def NN(self):
        del self._NN

    @property
    def bonds(self):
        """Return atom `Bonds` instance."""
        try:
            return Bonds([Bond(self, nn) for nn in self.NN])
        except (AttributeError, TypeError):
            return Bonds()

    @property
    def neighbors(self):
        """Neighbor atoms."""
        return self._neighbors

    @neighbors.setter
    def neighbors(self, value):
        if not isinstance(value, sknano.core.atoms.Atoms):
            raise TypeError('Expected an `Atoms` object.')
        self._neighbors = value

    @property
    def neighbor_distances(self):
        return self._neighbor_distances

    @neighbor_distances.setter
    def neighbor_distances(self, value):
        self._neighbor_distances = np.asarray(value)

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(CN=self.CN, NN=self.NN))
        return super_dict
