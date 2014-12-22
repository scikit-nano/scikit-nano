# -*- coding: utf-8 -*-
"""
==============================================================================
Nanotube bundle base class (:mod:`sknano.structures._nanotube_bundle`)
==============================================================================

.. currentmodule:: sknano.structures._nanotube_bundle

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from sknano.core.refdata import dVDW
from ._mixins import NanotubeBundleMixin

__all__ = ['NanotubeBundleBase']


class NanotubeBundleBase(NanotubeBundleMixin):
    """Nanotube bundle structure base class."""

    _bundle_geometries = ['square', 'rectangle', 'hexagon']

    def __init__(self, *Ch, nx=1, ny=1, vdw_spacing=dVDW,
                 bundle_packing=None, bundle_geometry=None, **kwargs):

        super(NanotubeBundleBase, self).__init__(*Ch, **kwargs)

        self.nx = nx
        self.ny = ny
        self.vdw_spacing = vdw_spacing

        self.bundle_geometry = bundle_geometry
        self.bundle_packing = bundle_packing

        self.generate_bundle_coords()
