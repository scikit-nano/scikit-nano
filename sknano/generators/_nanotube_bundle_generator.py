# -*- coding: utf-8 -*-
"""
===============================================================================
Nanotube bundle generator (:mod:`sknano.generators._nanotube_bundle_generator`)
===============================================================================

.. currentmodule:: sknano.generators._nanotube_bundle_generator

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import copy

__all__ = ['NanotubeBundleGeneratorBase']


class NanotubeBundleGeneratorBase:
    """Base class for generating nanotube bundles."""
    def __init__(self, autogen=True, **kwargs):

        super().__init__(autogen=False, **kwargs)

        if autogen:
            self.generate(generate_bundle_from_bundle_coords=True)

    def generate(self, generate_bundle=True,
                 generate_bundle_from_bundle_coords=False):
        """Generate structure data."""
        if generate_bundle and (generate_bundle_from_bundle_coords or
                                self.bundle_geometry is not None):
            super().generate()
            self.generate_bundle_from_bundle_coords()
        elif generate_bundle:
            self.generate_bundle()
            super().generate()
        else:
            super().generate()

    def generate_bundle(self):
        self.structure_data.clear()
        self.crystal_cell.scaling_matrix = [self.nx, self.ny, 1]

    def generate_bundle_from_bundle_coords(self):
        atomsobj0 = copy.deepcopy(self.atoms)
        atomsobj0.center_centroid()
        self.structure_data.clear()
        for mol_id, dr in enumerate(self.bundle_coords, start=1):
            atomsobj = copy.deepcopy(atomsobj0)
            atomsobj.translate(dr)
            [setattr(atom, 'mol', mol_id) for atom in atomsobj]
            self.atoms.extend(atomsobj)
            self.bundle_list.append(atomsobj)
