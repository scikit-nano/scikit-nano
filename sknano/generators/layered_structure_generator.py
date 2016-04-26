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
from operator import itemgetter
import importlib
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
    def parse_config(self):
        """Parse config file."""
        parser = self.parser
        parser.read(self.cfgfile)
        generator_module = 'sknano.generators'
        # generators = self.generators
        settings = self.settings

        fnames = []
        structures = []
        layers = []

        for section in parser.sections():
            if section == 'settings':
                [settings.update({option: value}) for option, value
                 in zip(parser[section].keys(), parser[section].values())]
                if self.verbose:
                    print(settings)
                continue
            else:
                layers.append(int(parser[section]['layer']))

            parameters = parser[section]['parameters']
            fname = '{}({})'.format(section[:-len('Generator')], parameters)
            fnames.append(fname.replace(' ', ''))

            call_sig = \
                self.call_signature.parseString(parameters, parseAll=True)[0]
            try:
                args, kwargs = call_sig
            except ValueError:
                args, kwargs = tuple(), call_sig[0]

            structure = getattr(importlib.import_module(generator_module),
                                section)(*args, **kwargs)
            structures.append(structure)
        self.fnames = [fname for (i, fname) in
                       sorted(zip(layers, fnames), key=itemgetter(0))]
        self.structures = [structure for (i, structure) in
                           sorted(zip(layers, structures), key=itemgetter(0))]

    def generate(self, finalize=True):
        """Generate structure data."""
        structures = self.structures
        settings = self.settings
        [structure.center_centroid() for structure in structures]
        for layer, structure in enumerate(structures[1:], start=1):
            # lattice_region = structure.lattice_region
            dy = -structure.bounding_box.ymin + \
                structures[layer-1].bounding_box.ymax - \
                float(settings.get('overlap', 0.0))
            structure.translate([0, dy, 0])

        atoms = StructureAtoms()
        [atoms.extend(structure.atoms) for structure in structures]
        atoms.center_centroid()
        atoms.assign_unique_ids()
        atoms.assign_unique_types()
        selstr = settings.get('selection', None)
        print('Natoms: {}'.format(atoms.Natoms))
        if selstr is not None:
            print('selstr: {}'.format(selstr))
            atoms = atoms.select(selstr)
            print('Natoms: {}'.format(atoms.Natoms))
        self.atoms.extend(atoms)
        # self.structure.extend(
        #     atoms.select(settings.get('selection', 'all')))

    def generate_fname(self):
        """Generate file name."""
        return '_on_'.join(self.fnames)
