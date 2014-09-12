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

from ._bonds import Bonds
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
    _atomattrs = XAtom._atomattrs + \
        ['CN', 'NN', 'bonds', 'pyramidalization_angle', 'sigma_bond_angle',
         'poav', 'poma']

    def __init__(self, **kwargs):
        super(KDTAtom, self).__init__(**kwargs)

        self._CN = 0
        self._NN = NeighborAtoms()
        self._bonds = Bonds()
        self._pyramidalization_angle = None
        self._sigma_bond_angle = None
        self.poav = None
        self.poma = []
        #self.mean_poma = None

    def __str__(self):
        """Return a nice string representation of `KDTAtom`."""
        strrep = "Atom(element={element!r}, atomID={atomID!r}, " + \
            "moleculeID={moleculeID!r}, atomtype={atomtype!r}, " + \
            "q={q!r}, m={m!r}, x={x:.6f}, y={y:.6f}, z={z:.6f}, " + \
            "CN={CN!r}, NN={NN!s}, bonds={bonds!s})"

        parameters = dict(element=self.element, atomID=self.atomID,
                          moleculeID=self.moleculeID, atomtype=self.atomtype,
                          q=self.q, m=self.m, x=self.x, y=self.y, z=self.z,
                          CN=self.CN, NN=self.NN, bonds=self.bonds)

        return strrep.format(**parameters)

    def __repr__(self):
        """Return canonical string representation of `KDTAtom`."""
        #strrep = "Atom(element={element!r}, atomID={atomID!r}, " + \
        #    "moleculeID={moleculeID!r}, atomtype={atomtype!r}, " + \
        #    "q={q!r}, m={m!r}, x={x:.6f}, y={y:.6f}, z={z:.6f}, " + \
        #    "CN={CN!r}, NN={NN!r})"
        #parameters = dict(element=self.element, atomID=self.atomID,
        #                  moleculeID=self.moleculeID, atomtype=self.atomtype,
        #                  q=self.q, m=self.m, x=self.x, y=self.y, z=self.z,
        #                  CN=self.CN, NN=self.NN)
        #return strrep.format(**parameters)
        return super(KDTAtom, self).__repr__()

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

    @property
    def bonds(self):
        return self._bonds

    @bonds.setter
    def bonds(self, value):
        if not isinstance(value, Bonds):
            raise TypeError('Expected a `Bonds` object.')
        self._bonds = value

    @property
    def sigma_bond_angle(self):
        return self._sigma_bond_angle

    @sigma_bond_angle.setter
    def sigma_bond_angle(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._sigma_bond_angle = value

    @property
    def pyramidalization_angle(self):
        return self._pyramidalization_angle

    @pyramidalization_angle.setter
    def pyramidalization_angle(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._pyramidalization_angle = value

    #@property
    #def poav_misalignment_angle(self):
    #    return self._poav_misalignment_angle

    #@poav_misalignment_angle.setter
    #def poav_misalignment_angle(self, value):
    #    if not isinstance(value, ):
    #        raise TypeError('Expected a number')
    #    self._poav_misalignment_angle = value
