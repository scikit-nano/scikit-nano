# -*- coding: utf-8 -*-
"""
===============================================================================
Nanotube generator base (:mod:`sknano.generators._base_nanotube_generator`)
===============================================================================

.. currentmodule:: sknano.generators._base_nanotube_generator

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import copy

from ._base import GeneratorBase

__all__ = ['NanotubeGeneratorBase']


class NanotubeGeneratorBase(GeneratorBase):
    """Base class for generating nanotubes."""
    def __init__(self, *args, autogen=True, generate_bundle=True,
                 from_bundle_coords=True, from_scaling_matrix=False, **kwargs):

        super().__init__(*args, autogen=False, **kwargs)

        if autogen:
            self.generate(generate_bundle=generate_bundle,
                          from_bundle_coords=from_bundle_coords,
                          from_scaling_matrix=from_scaling_matrix)

    def generate(self, generate_bundle=True, from_bundle_coords=False,
                 from_scaling_matrix=False):
        """Generate structure data."""
        if generate_bundle:
            if from_bundle_coords or self.bundle_geometry is not None:
                super().generate()
                self.generate_bundle_from_bundle_coords()
            else:
                self.generate_bundle_from_scaling_matrix()
                super().generate()
        else:
            super().generate()

    def generate_bundle_from_bundle_coords(self):
        """Generate bundle atoms by replicating structure atoms to \
            bundle coordinates."""
        atomsobj0 = copy.deepcopy(self.atoms)
        atomsobj0.center_centroid()
        self.structure_data.clear()
        for mol_id, dr in enumerate(self.bundle_coords, start=1):
            atomsobj = copy.deepcopy(atomsobj0)
            atomsobj.translate(dr)
            [setattr(atom, 'mol', mol_id) for atom in atomsobj]
            self.atoms.extend(atomsobj)
            self.bundle_list.append(atomsobj)

    def generate_bundle_from_scaling_matrix(self):
        """Generate bundle atoms by setting the \
            :attr:`~sknano.core.crystallography.CrystalCell.scaling_matrix`
            attribute."""
        self.structure_data.clear()
        self.crystal_cell.scaling_matrix = [self.nx, self.ny, 1]
