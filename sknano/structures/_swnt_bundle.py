# -*- coding: utf-8 -*-
"""
==============================================================================
SWNT bundle structure class (:mod:`sknano.structures._swnt_bundle`)
==============================================================================

.. currentmodule:: sknano.structures._swnt_bundle

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._compute_funcs import compute_bundle_density
from ._nanotube_bundle import NanotubeBundle
from ._swnt import SWNT

__all__ = ['SWNTBundle']


class SWNTBundle(NanotubeBundle, SWNT):

    def __init__(self, **kwargs):

        super(SWNTBundle, self).__init__(**kwargs)

    @property
    def bundle_density(self):
        return compute_bundle_density(self.n, self.m, d_vdw=self.vdw_spacing,
                                      bond=self.bond, element1=self.element1,
                                      element2=self.element2)

    @property
    def bundle_mass(self):
        return self.Ntubes * self.tube_mass

    @property
    def Natoms_per_bundle(self):
        return self.Ntubes * self.Natoms_per_tube

    @SWNT.Ntubes.getter
    def Ntubes(self):
        return len(self.bundle_coords)
