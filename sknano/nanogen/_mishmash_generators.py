# -*- coding: utf-8 -*-
"""
=================================================================
Mishmash generators (:mod:`sknano.nanogen._mishmash_generators`)
=================================================================

.. currentmodule:: sknano.nanogen._mishmash_generators

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#import copy
#import itertools
#import sys

#import numpy as np

from ..chemistry import Atoms

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

        See :py:meth:`~sknano.nanogen.StructureGenerator.save_data` method
        for documentation.

        """
        if fname is None:
            fname = 'mishmash_structure'
