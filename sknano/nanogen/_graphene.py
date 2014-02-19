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
import sys

import numpy as np

from pkshared.tools.refdata import CCbond

from ..chemistry import Atom, Atoms

__all__ = ['Graphene']


class Graphene(object):
    """Class for generating interactive Graphene objects."""

    def __init__(self, n=None, m=None, width=None, length=None, edge=None,
                 element1='C', element2='C', bond=CCbond, nlayers=1,
                 layer_spacing=3.35, stacking_order='AB',
                 autogen=True, with_units=False, verbose=False):

        self._n = n
        self._m = m

        self.element1 = element1
        self.element2 = element2

        self.Lx = width
        self.Ly = length
        self.edge = edge
        self.bond = bond
        self.verbose = verbose

        self.lx = 0.
        self.ly = 0.

        self.Nx = 0
        self.Ny = 0

        self.nlayers = nlayers
        self.layer_spacing = layer_spacing
        self.stacking_order = stacking_order

        self.layer_shift = np.zeros(3)

        if nlayers > 1 and stacking_order == 'AB':
            if edge in ('AC', 'armchair'):
                self.layer_shift[1] = self.bond
            elif edge in ('ZZ', 'zigzag'):
                self.layer_shift[0] = self.bond
            else:
                print('unrecognized edge parameter: {}'.format(edge))
                sys.exit(1)

        self.atom1 = Atom(element1)
        self.atom2 = Atom(element2)
        self.atom3 = Atom(element1)
        self.atom4 = Atom(element2)

        self._Natoms = 0
        self._Natoms_per_layer = None

        self.atoms = Atoms(atoms=[self.atom1,
                                  self.atom2,
                                  self.atom3,
                                  self.atom4])
        self.structure_atoms = None

        self._a1 = np.zeros(2, dtype=float)
        self._a2 = np.zeros(2, dtype=float)

        self._a1[0] = self._a2[0] = 1.5 * CCbond
        self._a1[1] = np.sqrt(3) / 2 * CCbond
        self._a2[1] = -self._a1[1]
