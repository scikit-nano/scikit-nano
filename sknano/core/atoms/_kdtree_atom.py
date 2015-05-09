# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class for KDTree analysis (:mod:`sknano.core.atoms._kdtree_atom`)
===============================================================================

An `XAtom` class for structure analysis.

.. currentmodule:: sknano.core.atoms._kdtree_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import sknano.core.atoms
from ._bond import Bond
from ._bonds import Bonds
from ._extended_atom import XAtom

__all__ = ['KDTAtom']


class KDTAtom(XAtom):
    """An `Atom` class for KDTree analysis."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['NN', 'bonds'])
        return attrs

    @XAtom.CN.getter
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

    @property
    def bonds(self):
        """Return atom `Bonds` instance."""
        try:
            return Bonds(bonds=[Bond(self, nn) for nn in self.NN])
        except (AttributeError, TypeError):
            return Bonds()

    #def todict(self):
    #    return dict(element=self.element, id=self.id,
    #                mol=self.mol, type=self.type,
    #                q=self.q, mass=self.mass,
    #                x=self.x, y=self.y, z=self.z,
    #                vx=self.vx, vy=self.vy, vz=self.vz,
    #                fx=self.fx, fy=self.fy, fz=self.fz,
    #                nx=self.nx, ny=self.ny, nz=self.nz,
    #                pe=self.pe, ke=self.ke, etotal=self.etotal,
    #                CN=self.CN, NN=self.NN, bonds=self.bonds)
