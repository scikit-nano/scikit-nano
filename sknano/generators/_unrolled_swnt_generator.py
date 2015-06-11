# -*- coding: utf-8 -*-
"""
===============================================================================
Unrolled SWNT generator (:mod:`sknano.generators._unrolled_swnt_generator`)
===============================================================================

.. currentmodule:: sknano.generators._unrolled_swnt_generator

.. todo::

   Add methods to perform fractional translation and cartesian translation
   before structure generation.

.. todo::

   Handle different units in output coordinates.

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
from builtins import range
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core import pluralize
from sknano.core.math import Vector
from sknano.structures import UnrolledSWNT
#from sknano.utils.geometric_shapes import Cuboid
from ._base import Atom, Atoms, GeneratorBase

__all__ = ['UnrolledSWNTGenerator']


class UnrolledSWNTGenerator(UnrolledSWNT, GeneratorBase):
    """Class for generating unrolled nanotube structures.

    .. versionadded:: 0.2.23

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions
    basis : {:class:`python:list`}, optional
        List of :class:`python:str`\ s of element symbols or atomic number
        of the two atom basis (default: ['C', 'C'])

        .. versionadded:: 0.3.10

    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.Atom` 1 and 2

        .. deprecated:: 0.3.10
           Use `basis` instead

    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    Lx, Ly, Lz : float, optional
        Length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.
    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.
    autogen : bool, optional
        if `True`, automatically call
        :meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :meth:`~NanotubeGenerator.generate_structure_data`.
    verbose : bool, optional
        if `True`, show verbose output

    Notes
    -----
    The `UnrolledSWNTGenerator` class generates graphene using the
    nanotube unit cell defined by the chiral vector
    :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    If you want to generate graphene with an armchair or zigzag edge using
    `length` and `width` parameters, see the
    :class:`~sknano.generators.GrapheneGenerator` class.

    .. seealso:: :class:`~sknano.generators.GrapheneGenerator`


    Examples
    --------
    First, load the :class:`~sknano.generators.UnrolledSWNTGenerator`
    class.

    >>> from sknano.generators import UnrolledSWNTGenerator

    Now let's generate an unrolled :math:`\\mathbf{C}_{\\mathrm{h}} = (10, 5)`
    SWCNT unit cell.

    >>> flatswcnt = UnrolledSWNTGenerator(n=10, m=5)
    >>> flatswcnt.save()

    The rendered structure looks like:

    """

    def __init__(self, *Ch, autogen=True, **kwargs):

        super().__init__(*Ch, **kwargs)

        if autogen:
            self.generate_structure_data()

    def generate_structure_data(self):
        """Generate structure data."""

        self.structure_data.clear()
        for nx in range(self.nx):
            for nz in range(int(np.ceil(self.nz))):
                dr = Vector([nx * self.Ch, 0.0, nz * self.T])
                for uc_atom in self.unit_cell:
                    nt_atom = Atom(element=uc_atom.symbol)
                    nt_atom.r = uc_atom.r + dr
                    self.atoms.append(nt_atom)

    def save(self, fname=None, outpath=None, structure_format=None,
             center_CM=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            chirality = '{}{}'.format('{}'.format(self.n).zfill(2),
                                      '{}'.format(self.m).zfill(2))

            nx = '{}' if self._assert_integer_nx else '{:.2f}'
            nx = ''.join((nx.format(self.nx), pluralize('cell', self.nx)))

            nz = '{}' if self._assert_integer_nz else '{:.2f}'
            nz = ''.join((nz.format(self.nz), pluralize('cell', self.nz)))

            cells = 'x'.join((nx, nz))
            fname_wordlist = (chirality, cells)

            fname = 'unrolled_' + '_'.join(fname_wordlist)

        if center_CM:
            self.atoms.center_CM()

        super().save(fname=fname, outpath=outpath,
                     structure_format=structure_format,
                     center_CM=False, **kwargs)
