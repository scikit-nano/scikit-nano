# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class for structure analysis (:mod:`sknano.core.atoms._structure_atom`)
===============================================================================

.. currentmodule:: sknano.core.atoms._structure_atom

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from ._cn_atom import CNAtom
from ._id_atom import IDAtom
from ._xyz_atom import XYZAtom
from ._charged_atom import ChargedAtom
from ._image_atom import ImageAtom
from ._type_atom import TypeAtom
from ._velocity_atom import VelocityAtom
from ._kdtree_atom import KDTreeAtomMixin
from ._poav_atom import POAVAtomMixin
from ._neighbor_atom import NeighborAtomMixin

__all__ = ['StructureAtom']


class StructureAtom(NeighborAtomMixin, POAVAtomMixin, KDTreeAtomMixin, CNAtom,
                    VelocityAtom, ImageAtom, XYZAtom, ChargedAtom, TypeAtom,
                    IDAtom):
    """An `Atom` class for structure analysis.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    id : int, optional
        atom ID
    mol : int, optional
        molecule ID
    type : int, optional
        atom type
    x, y, z : float, optional
        :math:`x, y, z` components of `StructureAtom` position vector relative
        to origin.
    vx, vy, vz : float, optional
        :math:`v_x, v_y, v_z` components of `StructureAtom` velocity.

    """
    def __init__(self, *args, NN=None, **kwargs):
        super().__init__(*args, **kwargs)

        self._neighbors = None
        if NN is not None:
            self.NN = NN

        self._POAV1 = None
        self._POAV2 = None
        self._POAVR = None
        # self.fmtstr = super().fmtstr + ", NN={NN!r}"

    def __dir__(self):
        attrs = super().__dir__()
        attrs.append('NN')
        # attrs.extend(['NN', 'bonds'])
        return attrs

    @CNAtom.CN.getter
    def CN(self):
        """`StructureAtom` coordination number."""
        try:
            return self.NN.Natoms
        except AttributeError:
            return super().CN
