# -*- coding: utf-8 -*-
"""
===============================================================================
Structure Generator base classes (:mod:`sknano.generators._base`)
===============================================================================

.. currentmodule:: sknano.generators._base

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import copy
import os

try:
    from sknano.core.atoms import KDTAtom as Atom, KDTAtoms as Atoms
except ImportError:
    from sknano.core.atoms import XAtom as Atom, XAtoms as Atoms

from sknano.io import StructureWriter, default_structure_format, \
    supported_structure_formats

__all__ = ['Atom', 'Atoms', 'GeneratorMixin', 'STRUCTURE_GENERATORS']


STRUCTURE_GENERATORS = ('SWNTGenerator', 'SWNTBundleGenerator',
                        'MWNTGenerator', 'MWNTBundleGenerator',
                        'GrapheneGenerator', 'BilayerGrapheneGenerator')


class GeneratorMixin(object):
    """Base mixin class for generator classes"""

    @property
    def atoms(self):
        """Return structure :class:`~sknano.core.atoms.Atoms`."""
        return self._atoms

    @atoms.setter
    def atoms(self, value):
        if not isinstance(value, Atoms):
            raise TypeError('Expected an `Atoms` object.')
        self._atoms = value

    @atoms.deleter
    def atoms(self):
        del self._atoms

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

        if outpath is not None:
            self.fpath = os.path.join(outpath, fname)

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
