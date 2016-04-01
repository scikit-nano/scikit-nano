# -*- coding: utf-8 -*-
"""
===================================================================================
Layered structure generator (:mod:`sknano.generators._layered_structure_generator`)
===================================================================================

.. currentmodule:: sknano.generators._layered_structure_generator

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

__all__ = ['LayeredStructureGenerator']

from collections import OrderedDict
import importlib
# import numpy as np

from sknano.core import BaseClass, call_signature
from sknano.core.atoms import StructureAtoms
from sknano.core.structures import StructureBase
from ._base import GeneratorBase

import configparser


class LayeredStructureGenerator(GeneratorBase, StructureBase, BaseClass):
    """Class for generating structures.

    Parameters
    ----------
    cfgfile : :class:`~python:str`

    """

    # call_signature = Forward()
    # args = Group(Optional(delimitedList(expr_item) + COMMA))
    # parameters = Group(args + ZeroOrMore(kwarg_expr))
    call_signature = call_signature.copy()

    def __init__(self, cfgfile=None, **kwargs):
        super().__init__(autogen=False, **kwargs)
        self.cfgfile = cfgfile
        self.fmtstr = "{cfgfile!r}"

        self.parser = configparser.ConfigParser()
        self.config = OrderedDict()
        self.structures = []
        self.fnames = []

        if cfgfile is not None:
            self._parse_config()

        self.generate()

    def _parse_config(self):
        parser = self.parser
        parser.read(self.cfgfile)

        del self.structures[:]
        for section in parser.sections():
            if section == 'settings':
                [self.config.update({option: value}) for option, value
                 in zip(parser[section].keys(), parser[section].values())]
                if self.verbose:
                    print(self.config)
                continue

            # generator_module = parser[section]['generator_module']
            generator_module = 'sknano.generators'
            generator_class = parser[section]['generator_class']
            parameters = parser[section]['parameters']
            fname = '{}({})'.format(generator_class[:-len('Generator')],
                                    parameters)
            self.fnames.append(fname.replace(' ', ''))
            # self.config.update({section: {'generator_class': generator_class,
            #                               'parameters': parameters}})

            call_sig = \
                self.call_signature.parseString(parameters, parseAll=True)[0]
            try:
                args, kwargs = call_sig
            except ValueError:
                args, kwargs = tuple(), call_sig[0]

            generator = getattr(importlib.import_module(generator_module),
                                generator_class)(*args, **kwargs)
            self.structures.append(generator)

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
        if selstr is not None:
            atoms = atoms.select(selstr)
        self.structure.extend(atoms)
        # self.structure.extend(
        #     atoms.select(self.config.get('selection', 'all')))

    def generate_fname(self):
        fname = '_on_'.join([fname for fname in reversed(self.fnames)])
        return fname

    def save(self, fname=None, structure_format=None, **kwargs):
        if fname is None:
            kwargs['fname'] = self.config.get('fname', self.generate_fname())
        if structure_format is None:
            kwargs['structure_format'] = \
                self.config.get('structure_format', 'xyz')
        super().save(**kwargs)

    def todict(self):
        return dict(cfgfile=self.cfgfile)
