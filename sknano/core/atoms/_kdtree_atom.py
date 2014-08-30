# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class for KDTree analysis (:mod:`sknano.core.atoms._kdtree_atom`)
===============================================================================

An `XAtom` class for structure analysis.

.. currentmodule:: sknano.core.atoms._kdtree_atom

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numbers

from ._extended_atom import XAtom
from ._neighbor_atoms import NeighborAtoms

__all__ = ['KDTAtom']


class KDTAtom(XAtom):
    """An `Atom` class for KDTree analysis.

    Parameters
    ----------
    CN : int, optional
        `KDTAtom` coordination number.
    NN : sequence, optional
        List of nearest-neighbor `KDTAtom` objects instances

    """

    def __init__(self, CN=0, NN=None, **kwargs):
        super(KDTAtom, self).__init__(**kwargs)

        self._CN = CN
        if NN is None or not isinstance(NN, NeighborAtoms):
            NN = NeighborAtoms()
        self._NN = NN

    def __repr__(self):
        """Return string representation of `KDTAtom`."""
        reprstr = "Atom(element={element!r}, atomID={atomID!r}, " + \
            "moleculeID={moleculeID!r}, atomtype={atomtype!r}, " + \
            "q={q!r}, m={m!r}, x={x:.6f}, y={y:.6f}, z={z:.6f}, " + \
            "CN={CN!r}, NN={NN!r})"

        parameters = dict(element=self.element, atomID=self.atomID,
                          moleculeID=self.moleculeID, atomtype=self.atomtype,
                          q=self.q, m=self.m, x=self.x, y=self.y, z=self.z,
                          CN=self.CN, NN=self.NN)

        return reprstr.format(**parameters)

    @property
    def CN(self):
        """Return `KDTAtom` coordination number."""
        return self._CN

    @CN.setter
    def CN(self, value):
        """Set `KDTAtom` coordination number."""
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number.')
        self._CN = int(value)

    @property
    def NN(self):
        """Return `KDTAtom` `NeighborAtoms` object."""
        return self._NN

    @NN.setter
    def NN(self, value):
        """Set `KDTAtom` `NeighborAtoms`."""
        if not isinstance(value, NeighborAtoms):
            raise TypeError('Expected `NeighborAtoms`.')
        self._NN = value

    def compute_sigma_bond_angles(self):
        pass

    def compute_pyramidalization_angle(self):
        pass
