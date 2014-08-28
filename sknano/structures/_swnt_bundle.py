# -*- coding: utf-8 -*-
"""
==============================================================================
SWNT bundle structure class (:mod:`sknano.structures._swnt_bundle`)
==============================================================================

.. currentmodule:: sknano.structures._swnt_bundle

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ..core.math import Vector
from ..core.refdata import dVDW

from ._mixins import NanotubeBundleMixin
from ._swnt import SWNT

__all__ = ['SWNTBundle']


class SWNTBundle(NanotubeBundleMixin, SWNT):

    def __init__(self, nx=1, ny=1, Lx=None, Ly=None, vdw_spacing=dVDW,
                 bundle_packing=None, bundle_geometry=None, **kwargs):

        super(SWNTBundle, self).__init__(**kwargs)

        self._nx = int(nx)
        self._ny = int(ny)
        self._Lx = Lx
        self._Ly = Ly

        self._r1 = Vector()

        self._r2 = Vector()

        self._vdw_spacing = vdw_spacing
        self._bundle_packing = bundle_packing
        self._bundle_geometry = bundle_geometry

        self._Natoms_per_bundle = None

        self._bundle_mass = None
        self._bundle_density = None

        self.compute_bundle_params()

    def compute_bundle_params(self):
        super(SWNTBundle, self).compute_tube_params()
        super(SWNTBundle, self).compute_bundle_params()
