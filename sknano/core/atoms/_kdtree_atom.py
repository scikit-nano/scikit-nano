# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class for KDTree analysis (:mod:`sknano.core.atoms._kdtree_atom`)
===============================================================================

An `XAtom` sub-class for structure analysis.

.. currentmodule:: sknano.core.atoms._kdtree_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._extended_atom import XAtom
from ._mixins import NNAtomMixin

__all__ = ['KDTAtom']


class KDTAtom(XAtom, NNAtomMixin):
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
            return super().NN.Natoms
        except AttributeError:
            return super().CN

    # def todict(self):
    #    return dict(element=self.element, id=self.id,
    #                mol=self.mol, type=self.type,
    #                q=self.q, mass=self.mass,
    #                x=self.x, y=self.y, z=self.z,
    #                vx=self.vx, vy=self.vy, vz=self.vz,
    #                fx=self.fx, fy=self.fy, fz=self.fz,
    #                nx=self.nx, ny=self.ny, nz=self.nz,
    #                pe=self.pe, ke=self.ke, etotal=self.etotal,
    #                CN=self.CN, NN=self.NN, bonds=self.bonds)
