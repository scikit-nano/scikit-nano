# -*- coding: utf-8 -*-
"""
============================================================================
Fullerene structure generators (:mod:`sknano.nanogen._fullerene_generators`)
============================================================================

.. currentmodule:: sknano.nanogen._fullerene_generators

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#import copy
#import itertools
#import sys

#import numpy as np

from ..chemistry import Atoms
#from ..tools.refdata import CCbond

from ._fullerenes import Fullerene
from ._structure_generator import StructureGenerator

__all__ = ['FullereneGenerator']


class FullereneGenerator(Fullerene, StructureGenerator):
    u"""Class for generating fullerene structures.

    .. note::

       This class is under development and not functional! Attempting
       to use it raises a `RuntimeError`.

    Parameters
    ----------
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    autogen : bool, optional
        if `True`, automatically call
        :py:meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :py:meth:`~NanotubeGenerator.generate_structure_data`.
    verbose : bool, optional
        if `True`, show verbose output

    Raises
    ------
    `RuntimeError`
        If called. This class is under development and not yet functional.

    Examples
    --------
    First, load the :py:class:`~sknano.nanogen.FullereneGenerator` class.

    >>> from sknano.nanogen import FullereneGenerator
    >>> fg = FullereneGenerator(N=60)

    """

    def __init__(self, N=int, PG=None, autogen=True, verbose=False):

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

        super(FullereneGenerator, self).save_data(
            fname=fname, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            deg2rad=deg2rad, center_CM=center_CM)
