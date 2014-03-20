# -*- coding: utf-8 -*-
"""
===============================================================================
Fullerene lattice (:mod:`sknano.nanogen._fullerene_bravais_lattice_generators`)
===============================================================================

.. currentmodule:: sknano.nanogen._fullerene_bravais_lattice_generators

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#import copy
#import itertools
#import sys

#import numpy as np

from ..chemistry import Atoms
#from ..tools.refdata import CCbond

from ._fullerene_generators import FullereneGenerator
#from ._lattice_generators import BravaisLatticeGenerator

__all__ = ['FullereneBravaisLatticeGenerator']


class FullereneBravaisLatticeGenerator(FullereneGenerator):
    u"""Class for generating fullerene bravais lattice structures.

    .. note::

       This class is under development and not functional! Attempting
       to use it raises a `RuntimeError`.

    """

    def __init__(self, n=int, autogen=True, verbose=False):

        raise RuntimeError('This class is not yet implemented.')

        super(FullereneGenerator, self).__init__(with_units=False,
                                                 verbose=verbose)

        if autogen:
            self.generate_structure_data()

    def generate_structure_data(self):
        """Generate structure data."""
        self._structure_atoms = Atoms()

    def save_data(self, fname=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True):
        """Save structure data.

        See :py:meth:`~sknano.nanogen.StructureGenerator.save_data` method
        for documentation.

        """
        if fname is None:
            fname = 'C{}'.format(self.Natoms)

        super(FullereneBravaisLatticeGenerator, self).save_data(
            fname=fname, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            deg2rad=deg2rad, center_CM=center_CM)
