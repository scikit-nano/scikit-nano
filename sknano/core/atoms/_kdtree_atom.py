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

import sknano.core.atoms
from ._bonds import Bonds
from ._extended_atom import XAtom

__all__ = ['KDTAtom']


class KDTAtom(XAtom):
    """An `Atom` class for KDTree analysis.

    Parameters
    ----------
    CN : int, optional
        `KDTAtom` coordination number.
    NN : {sequence, `Atoms`}, optional
        List of nearest-neighbor `Atoms`
    bonds : {sequence, `Bonds`}, optional
        List of atom `Bond`\ s or `Bonds` instance

    """
    _atomattrs = XAtom._atomattrs + ['CN', 'NN', 'bonds']

    def __init__(self, CN=0, NN=None, bonds=None, **kwargs):
        super(KDTAtom, self).__init__(**kwargs)

        self.CN = CN

        if NN is None:
            NN = sknano.core.atoms.StructureAtoms()
        self.NN = NN

        if bonds is None:
            bonds = Bonds()
        self.bonds = bonds

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
        """Nearest-neighbor `Atoms`."""
        return self._NN

    @NN.setter
    def NN(self, value):
        """Set nearest-neighbor `Atoms`."""
        if not isinstance(value, sknano.core.atoms.Atoms):
            raise TypeError('Expected an `Atoms` object.')
        self._NN = value

    @property
    def bonds(self):
        """Return atom `Bonds` instance."""
        return self._bonds

    @bonds.setter
    def bonds(self, value):
        if not isinstance(value, Bonds):
            raise TypeError('Expected a `Bonds` object.')
        self._bonds = value
