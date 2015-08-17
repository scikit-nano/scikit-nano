# -*- coding: utf-8 -*-
"""
===============================================================================
SWNT structure generator (:mod:`sknano.generators._swnt_generator`)
===============================================================================

.. currentmodule:: sknano.generators._swnt_generator

.. todo::

   Add methods to perform fractional translation and cartesian translation
   before structure generation.

.. todo::

   Handle different units in output coordinates.

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core import pluralize
from sknano.core.math import Point
from sknano.structures import SWNT
from sknano.core.geometric_regions import Cuboid
from ._base import GeneratorBase

__all__ = ['SWNTGenerator']


class SWNTGenerator(GeneratorBase, SWNT):
    """Class for generating :class:`SWNT` structures.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of integers (i.e., *Ch = [(n, m)]) or
        2 integers (i.e., *Ch = [n, m]) specifying the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nz : :class:`python:int`, optional
        Number of repeat unit cells in the :math:`z` direction, along
        the *length* of the nanotube.
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
    Lz : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the `nz` value.

        .. versionadded:: 0.2.5

    tube_length : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the `nz` value.

        .. deprecated:: 0.2.5
           Use `Lz` instead

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    autogen : bool, optional
        if `True`, automatically generate structure data.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------
    First, load the :class:`~sknano.generators.SWNTGenerator` class.

    >>> from sknano.generators import SWNTGenerator

    Now let's generate a :math:`\\mathbf{C}_{\\mathrm{h}} = (10, 5)`
    :class:`SWNT` unit cell.

    >>> swnt = SWNTGenerator((10, 5))
    >>> swnt.save(fname='10,5_unit_cell.xyz')
    >>> # note that there are two other alternative, but equivalent
    >>> # means of passing arguments to SWNTGenerator constructor:
    >>> # SWNTGenerator(10, 5) and SWNTGenerator(n=10, m=5)

    Here's a nice ray traced rendering of the generated
    :math:`(10, 5)` :class:`SWNT` unit cell.

    .. image:: /images/10,5_unit_cell-01.png

    """
    def generate(self):
        """Generate structure data."""

        super().generate()

        if self.L0 is not None and self.fix_Lz:
            Lz_cutoff = 10 * self.L0 + 1
            region_bounds = Cuboid(pmin=Point([-np.inf, -np.inf, 0]),
                                   pmax=Point([np.inf, np.inf, Lz_cutoff]))
            self.atoms.clip_bounds(region_bounds)

    @classmethod
    def generate_fname(cls, n=None, m=None, nz=None, fix_Lz=False,
                       **kwargs):

        chirality = '{}{}'.format('{}'.format(n).zfill(2),
                                  '{}'.format(m).zfill(2))
        nz_fmtstr = '{:.2f}' if fix_Lz else '{:.0f}'
        nz = ''.join((nz_fmtstr.format(nz), pluralize('cell', nz)))
        fname_wordlist = (chirality, nz)
        fname = '_'.join(fname_wordlist)
        return fname

    def save(self, fname=None, outpath=None, structure_format=None,
             center_CM=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            fname = self.generate_fname(n=self.n, m=self.m, nz=self.nz,
                                        fix_Lz=self.fix_Lz)

        super().save(fname=fname, outpath=outpath,
                     structure_format=structure_format,
                     center_CM=center_CM, **kwargs)
