# -*- coding: utf-8 -*-
"""
===============================================================================
Graphene structure classes (:mod:`sknano.structures._graphene`)
===============================================================================

.. currentmodule:: sknano.structures._graphene

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#import itertools
import numbers

import numpy as np

#from sknano.core.atoms import Atom
from sknano.core.math import Vector
from sknano.core.refdata import dVDW  # , grams_per_Da
from ._base import StructureBase
from ._extras import edge_types

__all__ = ['GraphenePrimitiveCell', 'Graphene']


class GraphenePrimitiveCell(StructureBase):
    """Graphene primitive unit cell structure class.

    Parameters
    ----------
    bond : float, optional
        bond length between nearest-neighbor atoms in **Angstroms**.

    """
    def __init__(self, **kwargs):
        super(GraphenePrimitiveCell, self).__init__(**kwargs)
        self.a = np.sqrt(3) * self.bond

        self.a1 = Vector(nd=2)
        self.a2 = Vector(nd=2)

        self.a1.x = self.a2.x = np.sqrt(3) / 2 * self.a
        self.a1.y = 1 / 2 * self.a
        self.a2.y = -self.a1.y

        self.b1 = Vector(nd=2)
        self.b2 = Vector(nd=2)

        self.b1.x = self.b2.x = 1 / np.sqrt(3) * 2 * np.pi / self.a
        self.b1.y = 2 * np.pi / self.a
        self.b2.y = -self.b1.y


class Graphene(StructureBase):
    """Graphene structure class.

    Parameters
    ----------
    length : float, optional
        Length of graphene sheet in **nanometers**
    width : float, optional
        Width of graphene sheet in **nanometers**
    edge : {'AC', 'armchair', 'ZZ', 'zigzag'}, optional
        **A**\ rm\ **C**\ hair or **Z**\ ig\ **Z**\ ag edge along
        the `length` of the sheet.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2
    bond : float, optional
        bond length between nearest-neighbor atoms in **Angstroms**.
    nlayers : int, optional
        Number of graphene layers.
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers
    verbose : bool, optional
        verbose output

    Notes
    -----
    For now, the graphene structure is generated using a
    conventional unit cell, not the primitive unit cell.

    .. todo::

       Add notes on unit cell calculation.

    """

    def __init__(self, length=None, width=None, edge=None, nlayers=1,
                 layer_spacing=dVDW, layer_rotation_angles=None,
                 stacking_order='AB', deg2rad=True, **kwargs):

        super(Graphene, self).__init__(**kwargs)

        self.length = length
        self.width = width
        if edge in ('armchair', 'zigzag'):
            edge = edge_types[edge]
        elif edge not in ('AC', 'ZZ'):
            print('unrecognized edge parameter: {}\n'.format(edge) +
                  'choosing one at random...\n')
            edge = np.random.choice(['AC', 'ZZ'])
            print('the randomly chosen edge type is: {}'.format(edge))
        self.edge = edge

        self._Nx = 0
        self._Ny = 0

        self.layer_mass = None
        self.Natoms = 0
        self.Natoms_per_layer = 0

        self.nlayers = nlayers
        self.layer_spacing = layer_spacing

        if layer_rotation_angles is not None and deg2rad:
            if isinstance(layer_rotation_angles, numbers.Number):
                layer_rotation_angles = np.radians(layer_rotation_angles)
            elif isinstance(layer_rotation_angles, (list, np.ndarray)):
                layer_rotation_angles = \
                    np.radians(np.asarray(layer_rotation_angles)).tolist()

        self.layer_rotation_angles = layer_rotation_angles
        self.stacking_order = stacking_order

        self.layer_shift = Vector()

        if nlayers > 1 and stacking_order == 'AB':
            if edge == 'AC':
                self.layer_shift.y = self.bond
            else:
                self.layer_shift.x = self.bond

        self.cell = Vector()
        if edge == 'AC':
            # Set up the unit cell with the armchair edge aligned
            # along the `y`-axis.
            self.cell.x = np.sqrt(3) * self.bond
            self.cell.y = 3 * self.bond
        else:
            # Set up the unit cell with the zigzag edge aligned
            # along the `y`-axis.
            self.cell.x = 3 * self.bond
            self.cell.y = np.sqrt(3) * self.bond

        self.Nx = int(np.ceil(10 * self.width / self.cell.x))
        self.Ny = int(np.ceil(10 * self.length / self.cell.y))

    def __repr__(self):
        retstr = 'Graphene(length={!r}, width={!r}, edge={!r}, ' + \
            'element1={!r}, element2={!r}, bond={!r}, nlayers={!r}, ' + \
            'layer_spacing={!r}, layer_rotation_angles={!r}, ' + \
            'stacking_order={!r})'
        return retstr.format(self.length, self.width, self.edge, self.element1,
                             self.element2, self.bond, self.nlayers,
                             self.layer_spacing, self.layer_rotation_angles,
                             self.stacking_order)
