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

# from ._base import GeneratorBase

__all__ = ['NanotubeGeneratorBase', 'NanotubeBundleGeneratorMixin']


class NanotubeBundleGeneratorMixin:
    """Mixin class for generating nanotube bundles."""
    def generate(self, from_scaling_matrix=False):
        """Generate structure data."""
        if self.is_bundle:
            if from_scaling_matrix:
                # First, set the CrystalCell.scaling_matrix, then
                # generate StructureAtoms from the CrystalCell.BasisAtoms
                self.generate_bundle_from_scaling_matrix()
                super().generate()
            else:
                # First, generate the StructureAtoms from the
                # CrystalCell.BasisAtoms, then copy the StructureAtoms to
                # the bundle coordinates
                super().generate()
                self.generate_bundle_from_bundle_coords()
        else:
            super().generate()

    def generate_bundle_from_bundle_coords(self):
        """Generate bundle atoms by replicating structure atoms to \
            bundle coordinates."""
        atomsobj0 = copy.deepcopy(self.atoms)
        atomsobj0.center_centroid()
        self.structure.clear()
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
        self.structure.clear()
        self.crystal_cell.scaling_matrix = [self.nx, self.ny, 1]


class NanotubeGeneratorBase(NanotubeBundleGeneratorMixin):
    """Base class for generating nanotubes."""
    def __init__(self, *args, autogen=True, from_scaling_matrix=False,
                 **kwargs):
        super().__init__(*args, autogen=False, **kwargs)

        if autogen:
            self.generate(from_scaling_matrix=from_scaling_matrix)
