# -*- coding: utf-8 -*-
"""
==============================================================================
Nanotube bundle base class (:mod:`sknano.structures._nanotube_bundle`)
==============================================================================

.. currentmodule:: sknano.structures._nanotube_bundle

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.core.refdata import dVDW
from ._mixins import NanotubeBundleMixin

__all__ = ['NanotubeBundleBase']


class NanotubeBundleBase(NanotubeBundleMixin):
    """Nanotube bundle structure base class."""

    _bundle_geometries = ['square', 'rectangle', 'hexagon']

    def __init__(self, *args, nx=1, ny=1, vdw_spacing=dVDW,
                 bundle_packing=None, bundle_geometry=None, **kwargs):

        super().__init__(*args, **kwargs)

        self.nx = nx
        self.ny = ny
        self.vdw_spacing = vdw_spacing

        self.bundle_geometry = bundle_geometry
        self.bundle_packing = bundle_packing
        self.bundle_list = []
        self.generate_bundle_coords()

    def todict(self):
        attrdict = super().todict()
        attrdict.update(dict(nx=self.nx, ny=self.ny,
                             vdw_spacing=self.vdw_spacing,
                             bundle_packing=self.bundle_packing,
                             bundle_geometry=self.bundle_geometry))
        return attrdict
