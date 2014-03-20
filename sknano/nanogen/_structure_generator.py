# -*- coding: utf-8 -*-
"""
===================================================================
Structure generator (:mod:`sknano.nanogen._structure_generator`)
===================================================================

.. currentmodule:: sknano.nanogen._structure_generator

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import os

from ..chemistry import Atoms
from ..structure_io import DATAWriter, XYZWriter, default_structure_format, \
    supported_structure_formats
from ..tools import rotation_matrix

__all__ = ['StructureGenerator']


class StructureGenerator(object):
    """Base class for structure generator classes."""

    def __init__(self):
        self._fname = None
        self._fpath = None
        self._structure_atoms = Atoms()
        self._structure_format = None
        self._unit_cell = Atoms()

    @property
    def fname(self):
        """Structure file name."""
        return self._fname

    @property
    def fpath(self):
        """Structure file path."""
        return self._fpath

    @property
    def structure_atoms(self):
        """Return structure :class:`~sknano.chemistry.Atoms`."""
        return self._structure_atoms

    @property
    def structure_format(self):
        """Structure file format."""
        return self._structure_format

    @property
    def unit_cell(self):
        """Return unit cell :class:`~sknano.chemistry.Atoms`."""
        return self._unit_cell

    def save_data(self, fname=None, outpath=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True):
        """Save structure data.

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
        self._fname = fname

        if outpath is not None:
            self._fpath = os.path.join(outpath, fname)

        self._structure_format = structure_format

        if center_CM:
            self.structure_atoms.center_CM()

        if rotation_angle is not None:
            R_matrix = rotation_matrix(rotation_angle,
                                       rot_axis=rot_axis,
                                       deg2rad=deg2rad)
            self.structure_atoms.rotate(R_matrix)

        if structure_format == 'data':
            DATAWriter.write(fname=fname, outpath=outpath,
                             atoms=self.structure_atoms)
        else:
            XYZWriter.write(fname=fname, outpath=outpath,
                            atoms=self.structure_atoms)
