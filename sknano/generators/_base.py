# -*- coding: utf-8 -*-
"""
===============================================================================
Structure generator base module (:mod:`sknano.generators._base`)
===============================================================================

.. currentmodule:: sknano.generators._base

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext en'

import copy
import os

from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from sknano.io import StructureData, StructureWriter, \
    default_structure_format, supported_structure_formats

__all__ = ['Atom', 'Atoms', 'GeneratorBase', 'STRUCTURE_GENERATORS']


#: Tuple of structure generator classes.
STRUCTURE_GENERATORS = ('FullereneGenerator',
                        'GrapheneGenerator',
                        'BilayerGrapheneGenerator',
                        'UnrolledSWNTGenerator',
                        'SWNTGenerator',
                        'SWNTBundleGenerator',
                        'MWNTGenerator',
                        'MWNTBundleGenerator')


class GeneratorBase:
    """Base class for generator classes"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.structure_data = StructureData()

    def __getattr__(self, name):
        if name != '_structure_data':
            return getattr(self._structure_data, name)

    @property
    def structure_data(self):
        """Return :class:`~sknano.io.StructureData` instance."""
        return self._structure_data

    @structure_data.setter
    def structure_data(self, value):
        """Set :class:`~sknano.io.StructureData` instance."""
        if not isinstance(value, StructureData):
            raise TypeError('Expected a `StructureData` object.')
        self._structure_data = value

    @structure_data.deleter
    def structure_data(self):
        del self._structure_data

    @property
    def structure(self):
        """Alias to :attr:`~GeneratorBase.structure_data`."""
        return self.structure_data

    @structure.setter
    def structure(self, value):
        """Alias to :attr:`~GeneratorBase.structure_data`."""
        self.structure_data = value

    @structure.deleter
    def structure(self):
        del self.structure_data

    def save_data(self, fname=None, outpath=None, structure_format=None,
                  center_CM=True, **kwargs):
        """Save structure data.

        .. todo::

           Use the unit cell to set the bounds on output.

        Parameters
        ----------
        fname : {None, str}, optional
            file name string
        outpath : str, optional
            Output path for structure data file.
        structure_format : {None, str}, optional
            chemical file format of saved structure data. Must be one of:

                - `xyz`
                - `data`

            If `None`, then guess based on `fname` file extension or
            default to `xyz` format.
        center_CM : bool, optional
            Center center-of-mass on origin.

        """
        if fname.endswith(supported_structure_formats) and \
                structure_format is None:
            for ext in supported_structure_formats:
                if fname.endswith(ext):
                    structure_format = ext
                    break
        elif structure_format is None or \
            structure_format not in supported_structure_formats or \
            (not fname.endswith(supported_structure_formats) and
             structure_format not in supported_structure_formats):
            structure_format = default_structure_format

        if not fname.endswith(structure_format):
            fname += '.' + structure_format
        self.fname = fname

        if outpath is not None:
            fpath = os.path.join(outpath, fname)
        else:
            fpath = os.path.join(os.getcwd(), fname)
        self.fpath = fpath

        #self._structure_format = structure_format

        atoms = copy.deepcopy(self.atoms)

        if center_CM:
            atoms.center_CM()

        if kwargs:
            atoms.rotate(**kwargs)

        StructureWriter.write(fname=fname, outpath=outpath,
                              structure_format=structure_format, atoms=atoms)
