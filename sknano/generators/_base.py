# -*- coding: utf-8 -*-
"""
===============================================================================
Structure generator base module (:mod:`sknano.generators._base`)
===============================================================================

.. currentmodule:: sknano.generators._base

"""
from __future__ import absolute_import, division, print_function
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


class GeneratorBase(object):
    """Base class for generator classes"""

    def __init__(self, **kwargs):
        self.structure_data = StructureData()
        super(GeneratorBase, self).__init__(**kwargs)

    @property
    def atoms(self):
        """Return structure :class:`~sknano.core.atoms.Atoms`."""
        return self.structure_data.atoms

    @atoms.setter
    def atoms(self, value):
        if not isinstance(value, Atoms):
            raise TypeError('Expected an `Atoms` object.')
        self.structure_data.atoms = value

    @atoms.deleter
    def atoms(self):
        del self.structure_data.atoms

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
                  rotation_angle=None, rot_axis=None, anchor_point=None,
                  deg2rad=True, center_CM=True, savecopy=True, **kwargs):
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
        rotation_angle : {None, float}, optional
            Angle of rotation
        rot_axis : {'x', 'y', 'z'}, optional
            Rotation axis
        anchor_point : array_like, optional
            Rotation axis origin
        deg2rad : bool, optional
            Convert `rotation_angle` from degrees to radians.
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

        if savecopy:
            atoms = copy.deepcopy(self.atoms)
        else:
            atoms = self.atoms

        if center_CM:
            atoms.center_CM()

        if rotation_angle is not None:
            atoms.rotate(angle=rotation_angle, rot_axis=rot_axis,
                         anchor_point=anchor_point, deg2rad=deg2rad)

        StructureWriter.write(fname=fname, outpath=outpath, atoms=atoms,
                              structure_format=structure_format)
