# -*- coding: utf-8 -*-
"""
===============================================================================
Graphene structure class (:mod:`sknano.structures._graphene`)
===============================================================================

.. currentmodule:: sknano.structures._graphene

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import itertools
import numbers

import numpy as np

# from sknano.core.atoms import Atom
from sknano.core.crystallography import Crystal2DLattice, UnitCell
from sknano.core.math import Vector
from sknano.core.refdata import aCC, dVDW  # , grams_per_Da
from ._base import StructureBase

__all__ = ['GraphenePrimitiveCell', 'GrapheneConventionalCell', 'Graphene']


class GraphenePrimitiveCell(UnitCell):
    """Graphene primitive unit cell structure class.

    Parameters
    ----------

    """
    def __init__(self, bond=aCC, a=np.sqrt(3)*aCC, basis=['C', 'C'],
                 coords=[[0, 0, 0], [aCC, 0, 0]], cartesian=True):

        lattice = Crystal2DLattice(a=a, b=a, gamma=60)
        lattice.rotate(angle=-np.pi/6)

        super().__init__(lattice, basis, coords, cartesian)


class GrapheneConventionalCell(UnitCell):
    def __init__(self, bond=aCC, a=3*aCC, b=np.sqrt(3)*aCC, basis=4*['C'],
                 coords=[[0, 0, 0], [aCC, 0, 0],
                         [3/2*aCC, np.sqrt(3)/2*aCC, 0],
                         [5/2*aCC, np.sqrt(3)/2*aCC, 0]], cartesian=True):
        lattice = Crystal2DLattice.rectangular(a=a, b=b)
        super().__init__(lattice, basis, coords, cartesian)


class Graphene(StructureBase):
    """Graphene structure class.

    Parameters
    ----------
    armchair_edge_length : float, optional
        Length of armchair edge in **nanometers**

        .. versionadded:: 0.3.10

    zigzag_edge_length : float, optional
        Length of zigzag edge in **nanometers**

        .. versionadded:: 0.3.10

    length : float, optional
        Length of armchair edge in **nanometers**

        .. deprecated:: 0.3.10
           Use `armchair_edge_length` instead

    width : float, optional
        Width of graphene sheet in **nanometers**

        .. deprecated:: 0.3.10
           Use `zigzag_edge_length` instead

    edge : {'AC', 'armchair', 'ZZ', 'zigzag'}, optional
        **A**\ rm\ **C**\ hair or **Z**\ ig\ **Z**\ ag edge along
        the `length` of the sheet.

        .. deprecated:: 0.3.10
           No longer used!

    basis : {:class:`python:list`}, optional
        List of :class:`python:str`\ s of element symbols or atomic number
        of the two atom basis (default: ['C', 'C'])

        .. versionadded:: 0.3.10

    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.Atom` 1 and 2

        .. deprecated:: 0.3.10
           Use `basis` instead

    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    nlayers : int, optional
        Number of graphene layers (default: 1)
    layer_spacing : float, optional
        Distance between layers in **Angstroms** (default: 3.35).
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers.
    verbose : bool, optional
        verbose output

    Notes
    -----
    For now, the graphene structure is generated using a
    conventional unit cell, not the primitive unit cell.

    .. todo::

       Add notes on unit cell calculation.

    """

    def __init__(self, armchair_edge_length=None, zigzag_edge_length=None,
                 basis=['C', 'C'], bond=aCC, nlayers=1, layer_spacing=dVDW,
                 layer_rotation_angles=None, layer_rotation_increment=None,
                 stacking_order='AB', degrees=True, cartesian=True, **kwargs):

        if 'deg2rad' in kwargs:
            degrees = kwargs['deg2rad']
            del kwargs['deg2rad']

        if 'length' in kwargs:
            armchair_edge_length = kwargs['length']
            del kwargs['length']

        if 'width' in kwargs:
            zigzag_edge_length = kwargs['width']
            del kwargs['width']

        super().__init__(basis=basis, bond=bond, **kwargs)

        self.unit_cell = GrapheneConventionalCell(bond=self.bond,
                                                  basis=self.basis)

        self.armchair_edge_length = armchair_edge_length
        self.zigzag_edge_length = zigzag_edge_length

        self._Nx = 0
        self._Ny = 0

        self.layer_mass = None
        self.Natoms = 0
        self.Natoms_per_layer = 0

        self.nlayers = nlayers
        self.layer_spacing = layer_spacing

        if layer_rotation_increment is not None and \
                layer_rotation_angles is None:
            layer_rotation_angles = layer_rotation_increment * \
                np.arange(self.nlayers)
        else:
            layer_rotation_angles = np.zeros(self.nlayers)
            degrees = False

        if layer_rotation_angles is not None and degrees:
            if isinstance(layer_rotation_angles, numbers.Number):
                layer_rotation_angles = np.radians(layer_rotation_angles)
            elif isinstance(layer_rotation_angles, (list, np.ndarray)):
                layer_rotation_angles = \
                    np.radians(np.asarray(layer_rotation_angles)).tolist()

        self.layer_rotation_angles = layer_rotation_angles
        self.stacking_order = stacking_order

        self.layer_shift = Vector()

        if nlayers > 1 and stacking_order == 'AB':
            self.layer_shift.x = self.bond

        self.Nx = int(np.ceil(10 * self.armchair_edge_length / self.unit_cell.a))
        self.Ny = int(np.ceil(10 * self.zigzag_edge_length / self.unit_cell.b))
        self.fmtstr = 'armchair_edge_length={armchair_edge_length!r}, ' + \
            'zigzag_edge_length={zigzag_edge_length!r}, ' + \
            'bond={bond!r}, basis={basis!r}, nlayers={nlayers!r}, ' + \
            'layer_spacing={layer_spacing!r}, ' + \
            'layer_rotation_angles={layer_rotation_angles!r}, ' + \
            'stacking_order={stacking_order!r}'

    def todict(self):
        return dict(armchair_edge_length=self.armchair_edge_length,
                    zigzag_edge_length=self.zigzag_edge_length,
                    bond=self.bond, basis=self.basis,
                    nlayers=self.nlayers, layer_spacing=self.layer_spacing,
                    layer_rotation_angles=self.layer_rotation_angles,
                    stacking_order=self.stacking_order)
