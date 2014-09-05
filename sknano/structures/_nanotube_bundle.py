# -*- coding: utf-8 -*-
"""
==============================================================================
Nanotube bundle structure class (:mod:`sknano.structures._nanotube_bundle`)
==============================================================================

.. currentmodule:: sknano.structures._nanotube_bundle

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from sknano.core.math import Vector
from sknano.core.refdata import dVDW
from ._mixins import NanotubeBundleMixin

__all__ = ['NanotubeBundle']


class NanotubeBundle(NanotubeBundleMixin):

    def __init__(self, nx=1, ny=1, Lx=None, Ly=None, vdw_spacing=dVDW,
                 bundle_packing=None, bundle_geometry=None, **kwargs):

        super(NanotubeBundle, self).__init__(**kwargs)

        self.nx = nx
        self.ny = ny

        self.r1 = Vector()
        self.r2 = Vector()

        self.vdw_spacing = vdw_spacing
        self.bundle_geometry = bundle_geometry
        self._bundle_packing = bundle_packing

        self.bundle_coords = []
        self.compute_bundle_params()
