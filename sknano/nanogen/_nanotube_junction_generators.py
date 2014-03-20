# -*- coding: utf-8 -*-
"""
=============================================================================
Generate SWNT junctions (:mod:`sknano.nanogen._nanotube_junction_generators`)
=============================================================================

.. currentmodule:: sknano.nanogen._nanotube_junction_generators

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#import copy
#import itertools
#import sys

#import numpy as np

#from ..chemistry import Atom, Atoms
#from ..tools import plural_word_check
from ..tools.refdata import CCbond

from ._nanotube_generators import NanotubeGenerator

__all__ = ['NanotubeJunctionGenerator']


class NanotubeJunctionGenerator(NanotubeGenerator):
    u"""Class for generating nanotube junction structures.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nz : int, optional
        Number of repeat unit cells in the :math:`z` direction, along
        the *length* of the nanotube.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :py:class:`~sknano.chemistry.Atoms` 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    Lz : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the `nz` value.
    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.
    autogen : bool, optional
        if `True`, automatically call
        :py:meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :py:meth:`~NanotubeGenerator.generate_structure_data`.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    """
    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False,
                 autogen=True, verbose=False):

        raise RuntimeError('This class is not yet implemented.')

        super(NanotubeJunctionGenerator, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz,
            element1=element1, element2=element2,
            bond=bond, Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz,
            with_units=False, verbose=verbose)

        if autogen:
            super(NanotubeJunctionGenerator, self).generate_unit_cell()
            super(NanotubeJunctionGenerator, self).generate_structure_data()
