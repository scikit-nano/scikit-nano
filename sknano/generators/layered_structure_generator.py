# -*- coding: utf-8 -*-
"""
===================================================================================
Layered structure generator (:mod:`sknano.generators.layered_structure_generator`)
===================================================================================

.. currentmodule:: sknano.generators.layered_structure_generator

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

# from collections import OrderedDict
# import importlib
# import numpy as np

# from sknano.core import call_signature
from sknano.core.atoms import StructureAtoms
from .base import CompoundStructureGenerator

# import configparser

__all__ = ['LayeredStructureGenerator']


class LayeredStructureGenerator(CompoundStructureGenerator):
    """Class for generating structures.

    Parameters
    ----------
    cfgfile : :class:`~python:str`

    """
    def generate(self):
        structures = self.structures
        [structure.center_centroid() for structure in structures]
        for layer, structure in enumerate(structures[1:], start=1):
            # lattice_region = structure.lattice_region
            dy = -structure.bounding_box.ymin + \
                structures[layer-1].bounding_box.ymax - \
                float(self.config.get('overlap', 0.0))
            structure.translate([0, dy, 0])

        atoms = StructureAtoms()
        [atoms.extend(structure.atoms) for structure in structures]
        atoms.center_centroid()
        atoms.assign_unique_ids()
        atoms.assign_unique_types()
        selstr = self.config.get('selection', None)
        print('Natoms: {}'.format(atoms.Natoms))
        if selstr is not None:
            print('selstr: {}'.format(selstr))
            atoms = atoms.select(selstr)
            print('Natoms: {}'.format(atoms.Natoms))
        self.structure.extend(atoms)
        # self.structure.extend(
        #     atoms.select(self.config.get('selection', 'all')))

    def generate_fname(self):
        fname = '_on_'.join([fname for fname in reversed(self.fnames)])
        return fname
