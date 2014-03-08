# -*- coding: utf-8 -*-
"""
=================================================================
Graphene structure tools (:mod:`sknano.nanogen._graphene`)
=================================================================

.. currentmodule:: sknano.nanogen._graphene

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext'

#import itertools

try:
    from pint import UnitRegistry
    ureg = UnitRegistry()
    Qty = ureg.Quantity
except ImportError:
    Qty = None

import numpy as np

from ..chemistry import Atom
from ..tools import Vector2D, Vector3D
from ..tools.refdata import CCbond

edge_types = {'armchair': 'AC', 'zigzag': 'ZZ'}

__all__ = ['GraphenePrimitiveCell', 'Graphene']


class GraphenePrimitiveCell(object):
    """Class for generating interactive Graphene primitive unit cell.

    .. versionadded:: 0.2.23

    Parameters
    ----------
    bond : float, optional
        bond length between nearest-neighbor atoms in **Angstroms**.

    """
    def __init__(self, bond=CCbond, with_units=False, units=None):

        if bond is None:
            bond = CCbond

        if with_units and Qty is None:
            with_units = False
        self._with_units = with_units

        if with_units and units is None:
            units = 'angstroms'
        self._units = units

        if with_units and isinstance(bond, float):
            self._bond = Qty(bond, units)
        else:
            self._bond = bond

        try:
            self._a = np.sqrt(3) * self._bond.magnitude
        except AttributeError:
            self._a = np.sqrt(3) * self._bond

        self._a1 = Vector2D(with_units=with_units, units=units)
        self._a2 = Vector2D(with_units=with_units, units=units)

        self._a1.x = self._a2.x = np.sqrt(3) / 2 * self._a
        self._a1.y = 1 / 2 * self._a
        self._a2.y = -self._a1.y

        self._b1 = Vector2D(with_units=with_units, units='1/{}'.format(units))
        self._b2 = Vector2D(with_units=with_units, units='1/{}'.format(units))

        self._b1.x = self._b2.x = 1 / np.sqrt(3) * 2 * np.pi / self._a
        self._b1.y = 2 * np.pi / self._a
        self._b2.y = -self._b1.y

    @property
    def a(self):
        """Length of graphene unit cell vector."""
        return self._a

    @property
    def a1(self):
        """:math:`a_1` unit vector."""
        return self._a1

    @property
    def a2(self):
        """:math:`a_2` unit vector."""
        return self._a2

    @property
    def b1(self):
        """:math:`b_1` reciprocal lattice vector."""
        return self._b1

    @property
    def b2(self):
        """:math:`b_2` reciprocal lattice vector."""
        return self._b2


class Graphene(object):
    """Class for generating interactive Graphene objects.

    Parameters
    ----------
    width : float, optional
        Width of graphene sheet in **nanometers**
    length : float, optional
        Length of graphene sheet in **nanometers**
    edge : {'AC', 'armchair', 'ZZ', 'zigzag'}, optional
        **A**\ rm\ **C**\ hair or **Z**\ ig\ **Z**\ ag edge along
        the `length` of the sheet.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :py:class:`~sknano.chemistry.Atoms` 1 and 2
    bond : float, optional
        bond length between nearest-neighbor atoms in **Angstroms**.
    nlayers : int, optional
        Number of graphene layers.
    layer_spacing : float, optional
        Distance between layers in **Angstroms**.
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of graphene layers
    autogen : bool, optional
        automatically generate unit cell and full structure
    verbose : bool, optional
        verbose output

    Notes
    -----
    For now, the graphene structure is generated using a
    conventional unit cell, not the primitive unit cell.

    .. todo::

       Add notes on unit cell calculation.

    """

    def __init__(self, width=None, length=None, edge=None,
                 element1='C', element2='C', bond=CCbond, nlayers=1,
                 layer_spacing=3.35, stacking_order='AB', with_units=False,
                 units=None, verbose=False):

        # add each parameter in the order I want them to appear in
        # verbose output mode
        self._element1 = element1
        self._element2 = element2

        self._width = width
        self._length = length
        if edge in ('armchair', 'zigzag'):
            edge = edge_types[edge]
        elif edge not in ('AC', 'ZZ'):
            print('unrecognized edge parameter: {}\n'.format(edge) +
                  'choosing one at random...\n')
            edge = np.random.choice(['AC', 'ZZ'])
            print('the randomly chosen edge type is: {}'.format(edge))
        self._edge = edge

        if bond is None:
            bond = CCbond

        if with_units and Qty is None:
            with_units = False
        self._with_units = with_units

        if with_units and units is None:
            units = 'angstroms'
        self._units = units

        if with_units and isinstance(bond, float):
            self._bond = Qty(bond, units)
        else:
            self._bond = bond

        self._verbose = verbose

        self._cell = Vector2D(with_units=with_units, units=units)

        self._Nx = 0
        self._Ny = 0

        self._layer_mass = None
        self._Natoms = 0
        self._Natoms_per_layer = None

        self._nlayers = nlayers
        self._layer_spacing = layer_spacing
        self._stacking_order = stacking_order

        self._layer_shift = Vector3D(with_units=with_units, units=units)

        if nlayers > 1 and stacking_order == 'AB':
            if edge == 'AC':
                self._layer_shift.y = self._bond
            else:
                self._layer_shift.x = self._bond

        if edge == 'AC':
            # Set up the unit cell with the armchair edge aligned
            # along the `y`-axis.
            self._cell.x = np.sqrt(3) * bond
            self._cell.y = 3 * bond
        else:
            # Set up the unit cell with the zigzag edge aligned
            # along the `y`-axis.
            self._cell.x = 3 * bond
            self._cell.y = np.sqrt(3) * bond

        self._Nx = int(np.ceil(10 * self._width / self._cell.x))
        self._Ny = int(np.ceil(10 * self._length / self._cell.y))

    def compute_layer_params(self):
        pass

    @property
    def Natoms(self):
        """Number of :py:class:`~sknano.chemistry.Atoms` in **unit cell**"""
        return self._Natoms

    @property
    def Natoms_per_layer(self):
        """Number of :py:class:`~sknano.chemistry.Atoms` in **layer**."""
        return self._Natoms_per_layer

    @classmethod
    def compute_Natoms_per_layer(cls, n=None, m=None, element1=None,
                                 element2=None, **kwargs):
        return 0

    @property
    def layer_mass(self):
        """Graphene layer mass in grams."""
        return self._layer_mass

    @classmethod
    def compute_layer_mass(cls, n=None, m=None, width=None, length=None,
                           edge=None, element1=None, element2=None,
                           with_units=False, units='grams', magnitude=True):
        """Compute graphene layer mass in **grams**.

        Parameters
        ----------
        n, m : int, optional
        width, length : float, optional
        edge : {'AC', 'ZZ'}
        element1, element2 : str, optional
        with_units : bool, optional
        units : str, optional
        magnitude : bool, optional

        Returns
        -------
        mass : float

        """
        Natoms_per_layer = Graphene.compute_Natoms_per_layer(
            n=n, m=m, element1=element1, element2=element2)

        if element1 is None:
            element1 = 'C'
        if element2 is None:
            element2 = 'C'

        atom1 = Atom(element1)
        atom2 = Atom(element2)

        mass = Natoms_per_layer * (atom1.m + atom2.m) / 2

        return mass

    @property
    def nlayers(self):
        return self._nlayers

    @property
    def layer_spacing(self):
        return self._layer_spacing
