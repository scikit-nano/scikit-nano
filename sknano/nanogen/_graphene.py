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
import warnings
warnings.filterwarnings('ignore')  # to suppress the Pint UnicodeWarning

try:
    from pint import UnitRegistry
    ureg = UnitRegistry()
    Qty = ureg.Quantity
except ImportError:
    Qty = None

import numpy as np

from ..chemistry import Atom
from ..tools.refdata import CCbond

edge_types = {'armchair': 'AC', 'zigzag': 'ZZ'}

__all__ = ['Graphene']


class Graphene(object):
    """Class for generating interactive Graphene objects.

    Parameters
    ----------
    width : float
        Width of graphene sheet in **nanometers**
    length : float
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

    def __init__(self, n=None, m=None, width=None, length=None,
                 edge=None, element1='C', element2='C', bond=CCbond,
                 nlayers=1, layer_spacing=3.35, stacking_order='AB',
                 with_units=False, verbose=False):

        self._n = n
        self._m = m

        self._element1 = element1
        self._element2 = element2

        self._Lx = width
        self._Ly = length
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

        if with_units and isinstance(bond, float):
            self._bond = Qty(bond, 'angstroms')
        else:
            self._bond = bond

        self._verbose = verbose

        try:
            self._a = np.sqrt(3) * self._bond.magnitude
        except AttributeError:
            self._a = np.sqrt(3) * self._bond

        self._a1 = np.zeros(2, dtype=float)
        self._a2 = np.zeros(2, dtype=float)

        self._a1[0] = self._a2[0] = np.sqrt(3) / 2 * self._a
        self._a1[1] = 1 / 2 * self._a
        self._a2[1] = -self._a1[1]

        self._b1 = np.zeros(2, dtype=float)
        self._b2 = np.zeros(2, dtype=float)

        self._b1[0] = self._b2[0] = \
            1 / np.sqrt(3) * 2 * np.pi / self._a
        self._b1[1] = 2 * np.pi / self._a
        self._b2[1] = -self._b1[1]

        if with_units:
            self._a1 = Qty(self._a1, 'angstrom')
            self._a2 = Qty(self._a2, 'angstrom')

            self._b1 = Qty(self._b1, '1/angstrom')
            self._b2 = Qty(self._b2, '1/angstrom')

        self._lx = 0.
        self._ly = 0.

        self._Nx = 0
        self._Ny = 0

        self._nlayers = nlayers
        self._layer_spacing = layer_spacing
        self._stacking_order = stacking_order

        self._layer_shift = np.zeros(3)

        if nlayers > 1 and stacking_order == 'AB':
            if edge == 'AC':
                self._layer_shift[1] = self._bond
            else:
                self._layer_shift[0] = self._bond

        self._Natoms = 0
        self._Natoms_per_layer = None

        if edge == 'AC':
            # Set up the unit cell with the armchair edge aligned
            # along the `y`-axis.
            self._lx = np.sqrt(3) * bond
            self._ly = 3 * bond
        else:
            # Set up the unit cell with the zigzag edge aligned
            # along the `y`-axis.
            self._lx = 3 * bond
            self._ly = np.sqrt(3) * bond

        self._Nx = int(np.ceil(10 * self._Lx / self._lx))
        self._Ny = int(np.ceil(10 * self._Ly / self._ly))

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
