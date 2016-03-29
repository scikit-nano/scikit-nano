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

# import numpy as np

# from sknano.core.crystallography import SuperCell
from ._base import GeneratorBase

__all__ = ['NanotubeBundleGeneratorBase']


class NanotubeBundleGeneratorBase(GeneratorBase):
    """Base class for generating nanotubes."""
    def __init__(self, *args, from_scaling_matrix=False, **kwargs):
        self.from_scaling_matrix = from_scaling_matrix
        super().__init__(*args, **kwargs)

    def generate(self, finalize=True):
        """Generate structure data."""
        super().generate(finalize=False)
        if self.is_bundle:
            if self.from_scaling_matrix:
                # First, set the CrystalCell.scaling_matrix, then
                # generate StructureAtoms from the CrystalCell.BasisAtoms
                print('generating bundle form scaling matrix')
                self._generate_bundle_from_scaling_matrix()
            else:
                # First, generate the StructureAtoms from the
                # CrystalCell.BasisAtoms, then copy the StructureAtoms to
                # the bundle coordinates
                self._generate_bundle_from_bundle_coords()
        if finalize:
            self.finalize()

    def _generate_bundle_from_bundle_coords(self):
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

    def _generate_bundle_from_scaling_matrix(self):
        """Generate bundle atoms by setting the \
            :attr:`~sknano.core.crystallography.CrystalCell.scaling_matrix`."""
        # self.structure.clear()
        # nz = int(np.ceil(self.nz))
        self.crystal_cell.scaling_matrix = [self.nx, self.ny, 1]
