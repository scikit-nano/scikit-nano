# -*- coding: utf-8 -*-
"""
===============================================================================
Structure generator base module (:mod:`sknano.generators._base`)
===============================================================================

.. currentmodule:: sknano.generators._base

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod

# import copy
import os

import numpy as np

# from sknano.core.atoms import MDAtom as Atom, MDAtoms as Atoms
from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from sknano.io import StructureWriterMixin, supported_structure_formats

__all__ = ['Atom', 'Atoms', 'GeneratorBase', 'BaseGenerator', 'GeneratorMixin',
           'BulkGeneratorBase', 'BaseBulkGenerator',
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


class GeneratorBase(StructureWriterMixin, metaclass=ABCMeta):
    """Base structure generator class."""
    def __init__(self, *args, autogen=True, **kwargs):

        super().__init__(*args, **kwargs)

        if autogen:
            self.generate()

    @abstractmethod
    def generate(self):
        """Generate structure data."""
        return NotImplementedError

    def finalize(self):
        """Finalize structure data by assigning unique ids and types to \
            structure atoms."""
        self.assign_unique_ids()
        self.assign_unique_types()

    def save(self, *args, **kwargs):
        """An alias for :meth:`~sknano.io.StructureWriterMixin.write`."""
        super().write(*args, **kwargs)

BaseGenerator = GeneratorBase


class GeneratorMixin:
    """Mixin class with concrete implementation of \
        :meth:`~GeneratorBase.generate` method."""
    def generate(self):
        """Concrete implementation of :meth:`~GeneratorBase.generate` \
            method."""
        self.structure.clear()
        for atom in self.crystal_cell:
            self.atoms.append(Atom(**atom.todict()))
        self.finalize()


class BulkGeneratorBase(GeneratorMixin, GeneratorBase):
    """Base class for the *bulk structure generator* classes."""
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

BaseBulkGenerator = BulkGeneratorBase
