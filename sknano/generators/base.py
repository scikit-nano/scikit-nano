# -*- coding: utf-8 -*-
"""
===============================================================================
Structure generator base module (:mod:`sknano.generators.base`)
===============================================================================

.. currentmodule:: sknano.generators.base

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod
from collections import OrderedDict
import importlib
import os

import configparser

import numpy as np
# import copy
from sknano.core import BaseClass, call_signature
from sknano.core.atoms import StructureAtom as Atom
from sknano.core.structures import StructureBase
from sknano.io import StructureWriterMixin, supported_structure_formats


__all__ = ['GeneratorBase', 'GeneratorMixin',
           'CrystalStructureGenerator', 'NanoStructureGenerator',
           'CompoundStructureGenerator', 'STRUCTURE_GENERATORS']


#: Tuple of structure generator classes.
STRUCTURE_GENERATORS = ('FullereneGenerator',
                        'GrapheneGenerator',
                        'PrimitiveCellGrapheneGenerator',
                        'ConventionalCellGrapheneGenerator',
                        'BilayerGrapheneGenerator',
                        'UnrolledSWNTGenerator',
                        'SWNTGenerator',
                        'MWNTGenerator',
                        'AlphaQuartzGenerator',
                        'DiamondStructureGenerator',
                        'FCCStructureGenerator',
                        'CaesiumChlorideStructureGenerator',
                        'RocksaltStructureGenerator',
                        'ZincblendeStructureGenerator',
                        'MoS2Generator')


class GeneratorBase(StructureWriterMixin, StructureBase, metaclass=ABCMeta):
    """Base structure generator class.

    Parameters
    ----------
    autogen : :class:`~python:bool`
    finalize : :class:`~python:bool`

    """
    def __init__(self, *args, autogen=True, finalize=True, **kwargs):
        super().__init__(*args, **kwargs)

        if autogen:
            self.generate(finalize=finalize)

    # def __getattr__(self, name):
    #     print('getting attribute: {}'.format(name))
    #     if name != 'atoms' and len(self.atoms) > 0:
    #         attr = getattr(self.atoms, name, None)
    #         if attr:
    #             return attr
    #     super().__getattr__(name)

    @property
    def Natoms(self):
        """N atoms."""
        try:
            Natoms = self._atoms.Natoms
            if Natoms:
                return Natoms
            raise AttributeError
        except AttributeError:
            return super().Natoms

    @property
    def mass(self):
        """Total mass of atoms."""
        try:
            mass = self._atoms.mass
            if mass:
                return mass
            raise AttributeError
        except AttributeError:
            return super().mass

    @abstractmethod
    def generate(self, finalize=True):
        """Generate structure data."""
        return NotImplementedError

    def finalize(self):
        """Finalize structure data by assigning unique ids and types to \
            structure atoms."""
        atoms = self._atoms

        atoms.assign_unique_ids()
        atoms.assign_unique_types()
        self.structure.translate(self.lattice_shift)
        # atoms.lattice = self.crystal_cell.lattice

    def save(self, *args, **kwargs):
        """An alias for :meth:`~sknano.io.StructureWriterMixin.write`."""
        super().write(*args, **kwargs)


class GeneratorMixin:
    """Mixin class with concrete implementation of \
        :meth:`~GeneratorBase.generate` method."""
    def generate(self, finalize=True):
        """Concrete implementation of :meth:`~GeneratorBase.generate` \
            method."""
        self.structure.clear()
        for atom in self.crystal_cell:
            self.atoms.append(Atom(**atom.todict()))
        if finalize:
            self.finalize()


class CrystalStructureGenerator(GeneratorMixin, GeneratorBase):
    """:class:`GeneratorBase` sub-class for \
        :class:`~sknano.core.structures.CrystalStructureBase`\ s"""
    def save(self, fname=None, scaling_matrix=None, **kwargs):
        """Save structure data."""
        if fname is None:
            fname = self.__class__.__name__[:-len('Generator')]
            if fname.endswith('CC'):
                fname = '_'.join((fname, '-'.join(set(self.basis.symbols))))
        if scaling_matrix is None:
            scaling_matrix = self.scaling_matrix.A
        elif isinstance(scaling_matrix, np.matrix):
            scaling_matrix = scaling_matrix.A

        if fname is not None and scaling_matrix is not None:
            ext = None
            if fname.endswith(supported_structure_formats):
                fname, ext = os.path.splitext(fname)

            if isinstance(scaling_matrix, np.ndarray):
                if scaling_matrix.ndim == 2:
                    if scaling_matrix.shape == (1, 3):
                        scaling_matrix = scaling_matrix.flatten()
                    elif scaling_matrix.shape == (3, 3) and \
                            np.all(scaling_matrix -
                                   np.diag(np.diag(scaling_matrix)) == 0):
                        scaling_matrix = scaling_matrix.diagonal()
                    else:
                        scaling_matrix = scaling_matrix.flatten()

            fname = '_'.join((fname, 'x'.join(map(str, scaling_matrix))))
            if ext is not None:
                fname = fname + ext
        super().save(fname=fname, **kwargs)


class NanoStructureGenerator(GeneratorBase):
    """:class:`GeneratorBase` sub-class for \
        :class:`~sknano.core.structures.NanoStructureBase` class."""
    pass


class CompoundStructureGenerator(GeneratorBase, BaseClass):
    """Base class for compound structures.

    Parameters
    ----------
    cfgfile : :class:`~python:str`

    """
    # call_signature = Forward()
    # args = Group(Optional(delimitedList(expr_item) + COMMA))
    # parameters = Group(args + ZeroOrMore(kwarg_expr))
    call_signature = call_signature.copy()

    def __init__(self, cfgfile=None, autogen=True, **kwargs):
        super().__init__(autogen=False, **kwargs)
        self.fmtstr = "{cfgfile!r}"
        self.cfgfile = cfgfile
        self.parser = configparser.ConfigParser()
        self.config = OrderedDict()
        self.structures = []
        self.fnames = []

        if cfgfile is not None:
            self._parse_config()

        if autogen:
            self.generate()

    def _parse_config(self):
        parser = self.parser
        parser.read(self.cfgfile)
        generator_module = 'sknano.generators'

        del self.structures[:]
        fnames = []
        structures = []
        for section in parser.sections():
            if section == 'settings':
                [self.config.update({option: value}) for option, value
                 in zip(parser[section].keys(), parser[section].values())]
                if self.verbose:
                    print(self.config)
                continue

            parameters = parser[section]['parameters']
            fname = '{}({})'.format(section[:-len('Generator')], parameters)
            fnames.append(fname.replace(' ', ''))

            call_sig = \
                self.call_signature.parseString(parameters, parseAll=True)[0]
            try:
                args, kwargs = call_sig
            except ValueError:
                args, kwargs = tuple(), call_sig[0]

            generator = getattr(importlib.import_module(generator_module),
                                section)(*args, **kwargs)
            structures.append(generator)
        self.fnames = list(reversed(fnames[:]))
        self.structures = list(reversed(structures[:]))

    def save(self, fname=None, structure_format=None, **kwargs):
        if fname is None:
            kwargs['fname'] = self.config.get('fname', self.generate_fname())
        if structure_format is None:
            kwargs['structure_format'] = \
                self.config.get('structure_format', 'xyz')
        super().save(**kwargs)

    def todict(self):
        return dict(cfgfile=self.cfgfile)
