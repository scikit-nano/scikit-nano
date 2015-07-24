# -*- coding: utf-8 -*-
"""
==============================================================================
Generator mixin classes (:mod:`sknano.generators._mixins`)
==============================================================================

.. currentmodule:: sknano.generators._mixins

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import copy

__all__ = ['CappedNanotubeGeneratorMixin',
           'NanotubeBundleGeneratorMixin']


class CappedNanotubeGeneratorMixin:
    """Mixin class for generating capped nanotubes."""

    def generate_endcaps(self):
        pass


class NanotubeBundleGeneratorMixin:
    """Mixin class for generating nanotube bundles."""

    def generate_bundle(self):
        self.bundle_list = []
        atomsobj0 = copy.deepcopy(self.atoms)
        atomsobj0.center_CM()
        self.structure_data.clear()
        for mol_id, dr in enumerate(self.bundle_coords, start=1):
            atomsobj = copy.deepcopy(atomsobj0)
            atomsobj.translate(dr)
            [setattr(atom, 'mol', mol_id) for atom in atomsobj]
            self.atoms.extend(atomsobj)
            self.bundle_list.append(atomsobj)
