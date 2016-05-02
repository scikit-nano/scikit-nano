# -*- coding: utf-8 -*-
"""
========================================================================================
Nanotube bundle generator classes (:mod:`sknano.generators.nanotube_bundle_generator`)
========================================================================================

.. currentmodule:: sknano.generators.nanotube_bundle_generator

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import copy

# import numpy as np

# from sknano.core.crystallography import SuperCell
# from sknano.core.math import Vector
from sknano.core import grouper
from .base import NanoStructureGenerator

__all__ = ['NanotubeBundleGeneratorBase']


class NanotubeBundleGeneratorBase(NanoStructureGenerator):
    """Base class for generating nanotubes."""
    def generate(self, finalize=True):
        """Generate structure data."""
        super().generate(finalize=False)
        if self.is_bundle:
            if self.bundle_geometry is not None or \
                    self.__class__.__qualname__ == 'MWNTGenerator' or \
                    self.fix_L:
                self._generate_bundle_from_bundle_coords()
                # for mol, dr in enumerate(self.bundle_coords, start=1):
                #     atoms = self.atoms.filtered(self.atoms.mols == mol)
                #     atomsobj.translate(dr)
                #     [setattr(atom, 'mol', mol_id) for atom in atomsobj]
                #     self.bundle_list.append(atomsobj)
            else:
                atoms = self.atoms
                Natoms_per_tube = self.unit_cell.basis.Natoms * self.n3
                # if atoms.Natoms % Natoms_per_tube == 0:
                for mol, atoms in enumerate(grouper(atoms, Natoms_per_tube),
                                            start=1):
                    [setattr(atom, 'mol', mol) for atom in atoms]

        if finalize:
            self.finalize()

    def _generate_bundle_from_bundle_coords(self):
        """Generate bundle atoms by replicating structure atoms to \
            bundle coordinates."""
        atomsobj0 = copy.deepcopy(self.atoms)
        # atomsobj0 = copy.deepcopy(self.atoms.filtered(self.atoms.mols == 1))
        # atomsobj0.center_centroid()
        # self.lattice_shift = Vector(p0=self.crystal_cell.basis.centroid,
        #                             p=atomsobj0.centroid)

        self.structure.clear()
        for mol, dr in enumerate(self.bundle_coords, start=1):
            atomsobj = copy.deepcopy(atomsobj0)
            atomsobj.translate(dr)
            [setattr(atom, 'mol', mol) for atom in atomsobj]
            self.atoms.extend(atomsobj)
            self.bundle_list.append(atomsobj)
