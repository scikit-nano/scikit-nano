# -*- coding: utf-8 -*-
"""
=================================================================
Mishmash generators (:mod:`sknano.nanogen._mishmash_generators`)
=================================================================

.. currentmodule:: sknano.nanogen._mishmash_generators

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

#import copy
#import itertools
#import sys

#from pint import UnitRegistry
#ureg = UnitRegistry()
#Qty = ureg.Quantity

#import numpy as np

from ..chemistry import Atoms
from ..structure_io import DATAWriter, XYZWriter, default_structure_format, \
    supported_structure_formats
from ..tools import rotation_matrix

from ._structure_generator import StructureGenerator

__all__ = ['MishmashGenerator']


class MishmashGenerator(StructureGenerator):
    u"""Class for generating a mishmash of nanostructures.

    .. note::

       This class is in development and is not functional and will
       explode if you try and use it.

    Parameters
    ----------
    autogen : bool, optional
        if `True`, automatically call
        :py:meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :py:meth:`~NanotubeGenerator.generate_structure_data`.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    """

    def __init__(self):
        super(MishmashGenerator, self).__init__()

    def generate_structure_data(self):
        """Generate structure data."""
        self._structure_atoms = Atoms()

    def save_data(self, fname=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True):
        """Save structure data.

        Parameters
        ----------
        fname : {None, str}, optional
            file name string
        structure_format : {None, str}, optional
            chemical file format of saved structure data. Must be one of:

                - xyz
                - data

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
        if (fname is None and structure_format not in
                supported_structure_formats) or \
                (fname is not None and not
                    fname.endswith(supported_structure_formats) and
                    structure_format not in supported_structure_formats):
            structure_format = default_structure_format

        if fname is None:
            fname = 'mishmash_structure'
            fname += '.' + structure_format
        else:
            if fname.endswith(supported_structure_formats) and \
                    structure_format is None:
                for ext in supported_structure_formats:
                    if fname.endswith(ext):
                        structure_format = ext
                        break
            else:
                if structure_format is None or \
                        structure_format not in supported_structure_formats:
                    structure_format = default_structure_format

        self._fname = fname

        if center_CM:
            self._structure_atoms.center_CM()

        if rotation_angle is not None:
            R_matrix = rotation_matrix(rotation_angle,
                                       rot_axis=rot_axis,
                                       deg2rad=deg2rad)
            self._structure_atoms.rotate(R_matrix)

        if structure_format == 'data':
            DATAWriter.write(fname=self._fname, atoms=self._structure_atoms)
        else:
            XYZWriter.write(fname=self._fname, atoms=self._structure_atoms)
