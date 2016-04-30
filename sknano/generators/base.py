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
import os

import configparser

import numpy as np
# import copy
from sknano.core import BaseClass, call_signature
from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from sknano.core.structures import StructureBase
from sknano.io import StructureReaderMixin, StructureWriterMixin, \
    supported_structure_formats


__all__ = ['GeneratorBase', 'GeneratorMixin',
           'CrystalStructureGenerator', 'NanoStructureGenerator',
           'CompoundStructureGenerator', 'DefectGeneratorBase',
           'STRUCTURE_GENERATORS']


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

        self._atoms = Atoms()

        if autogen:
            self.generate(finalize=finalize)

    def __getattr__(self, name):
        # print('getting attribute: {}'.format(name))
        # if name != 'atoms' and len(self.atoms) > 0:
        #     attr = getattr(self.atoms, name, None)
        #     if attr:
        #         return attr
        try:
            return getattr(self.atoms, name)
        except AttributeError:
            return super().__getattr__(name)

    # def __eq__(self, other):
    #     if not isinstance(other, self.__class__):
    #         return NotImplemented
    #     return (self is other or self._atoms == other._atoms)

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
        self.generators = OrderedDict()
        self.settings = OrderedDict()
        self.structures = []
        self.fnames = []

        if cfgfile is not None:
            self.parse_config()

        if autogen:
            self.generate()

    @abstractmethod
    def parse_config(self):
        """Parse config file."""
        return NotImplementedError

    def save(self, fname=None, structure_format=None, **kwargs):
        """Save structure data."""
        settings = self.settings
        if fname is None:
            kwargs['fname'] = settings.get('fname', self.generate_fname())
        if structure_format is None:
            kwargs['structure_format'] = \
                settings.get('structure_format', 'xyz')
        super().save(**kwargs)

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(cfgfile=self.cfgfile)


class DefectGeneratorBase(StructureReaderMixin, GeneratorBase, BaseClass):
    """Base class for generating structure data with defects.

    Parameters
    ----------

    """
    def __init__(self, structure_obj=None, autogen=True, **kwargs):
        super().__init__(autogen=False, **kwargs)
        self.fmtstr = "{structure!r}"
        if structure_obj is not None:
            if not isinstance(structure_obj, StructureBase):
                structure_obj = self.read(structure_obj).structure
            else:
                structure_obj = structure_obj.structure
        self.structure = structure_obj
        self.__update_attrs()

    def update_attrs(self):
        structure = self.structure
        if structure is not None:
            self._atoms = structure.atoms
            self._crystal_cell = structure.crystal_cell
            self._lattice_shift = structure.lattice_shift
            self._region = structure.region

    __update_attrs = update_attrs

    @property
    def atom_ids(self):
        return self.atoms.ids

    @property
    def atom_coords(self):
        return self.atoms.get_coords(asdict=True)

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(structure=self.structure)
