# -*- coding: utf-8 -*-
"""
===============================================================================
MWNT structure generator (:mod:`sknano.generators._mwnt_generator`)
===============================================================================

.. currentmodule:: sknano.generators._mwnt_generator

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

# import copy

# import numpy as np

# from sknano.core import pluralize
# from sknano.core.math import Vector
from sknano.structures import MWNT
# from sknano.core.geometric_regions import Cuboid
from ._base import GeneratorBase
from ._swnt_generator import SWNTGenerator

__all__ = ['MWNTGenerator']


class MWNTGenerator(GeneratorBase, MWNT):
    """Class for generating single, `MWNT`.

    .. versionchanged:: 0.2.20

       `MWNTGenerator` no longer generates MWNT *bundles*, only *single*
       MWNTs. To generate bundled MWNT structure data, use the
       `MWNTBundleGenerator` class.

    .. versionadded:: 0.2.8

    Parameters
    ----------
    Ch_list : :class:`python:list`, optional
        (:attr:`~SWNT.n`, :attr:`~SWNT.m`) for each `SWNT` wall in `MWNT`.
    Nwalls : int, optional
        Number of `SWNT` walls in `MWNT`.
    Lz : float, optional
        `MWNT` length in **nanometers**.
    min_wall_diameter : float, optional
        Minimum `MWNT` wall diameter, in units of **Angstroms**.
    max_wall_diameter : float, optional
        Maximum `MWNT` wall diameter, in units of **Angstroms**.
    max_walls : int, optional
        Maximum number of `MWNT` walls.
    chiral_types : {None, 'armchair', 'zigzag', 'achiral', 'chiral'}, optional
        If `None`, the :attr:`~SWNT.chiral_type` of each `MWNT` walls
        will be random and determined by the set of randomly selected
        chiral indices (:attr:`~SWNT.n`, :attr:`~SWNT.m`).
    wall_spacing : float, optional
        Inter-wall spacing in units of **Angstroms**.
        Default value is the van der Waals interaction distance of 3.4
        Angstroms.
    autogen : bool, optional
        if `True`, automatically call
        :meth:`~MWNTGenerator.generate`.
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
        nearest neighbor atoms, in units of **Angstroms**.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.generators import MWNTGenerator
    >>> MWNTGenerator(Nwalls=5, min_wall_diameter=10, Lz=5).save()

    The above command generated a 5 wall,  10 **nanometer** long `MWNT`.
    The only constraints on the `MWNT` wall chiralities were the diameter
    constraints imposed by the `min_wall_diameter` parameter,
    which set the minimum wall diameter to 10 Angstroms,
    as well as the minimum wall-to-wall separation, which
    defaults to the van der Waals distance of 3.4 Angstroms.

    This `MWNT` chirality may be written as:

    :math:`\\mathbf{C}_{\\mathrm{h}} = (8,7)@(17,8)@(9,24)@(27,18)@(22,32)`

    Here's a colorful rendering of the generated `MWNT` structure:

    .. image:: /images/5wall_mwnt_(8,7)@(17,8)@(9,24)@(27,18)@(22,32)-04.png

    """
    def generate(self):
        """Generate structure data.

        .. todo::

           Load the diameter and chirality data from file instead of
           generating it every time.

        """

        self.structure_data.clear()
        for swnt in self.walls:
            self.atoms.extend(SWNTGenerator(**swnt.todict()).atoms)

    @classmethod
    def generate_fname(cls, Ch_list=None, Nwalls=None, **kwargs):
        Nwalls = '{}wall_mwnt'.format(len(Ch_list))
        chiralities = '@'.join([str(Ch).replace(' ', '') for
                                Ch in Ch_list])

        fname = '_'.join((Nwalls, chiralities))
        return fname

    def save(self, fname=None, outpath=None, structure_format=None,
             center_centroid=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            fname = self.generate_fname(Ch_list=self.Ch_list)

        super().save(fname=fname, outpath=outpath,
                     structure_format=structure_format,
                     center_centroid=center_centroid, **kwargs)
