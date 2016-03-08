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

from sknano.core import pluralize
# from sknano.core.math import Vector
from sknano.core.structures import MWNT
# from sknano.core.geometric_regions import Cuboid
from ._base import GeneratorBase
from ._nanotube_generator_base import NanotubeGeneratorBase
from ._swnt_generator import SWNTGenerator

__all__ = ['MWNTGeneratorMixin', 'MWNTGenerator']


class MWNTGeneratorMixin:
    """`MWNTGenerator` mixin."""
    def generate(self, **kwargs):
        """Generate structure data.

        .. todo::

           Load the diameter and chirality data from file instead of
           generating it every time.

        """
        self.structure.clear()
        for swnt in self.walls:
            self.atoms.extend(
                SWNTGenerator(fix_Lz=True, **swnt.todict()).atoms)
        self.finalize()


class MWNTGenerator(NanotubeGeneratorBase, MWNTGeneratorMixin,
                    GeneratorBase, MWNT):
    """Class for generating multi-walled nanotubes (MWNT).

    .. versionchanged:: 0.3.22

       `MWNTGenerator` now generates both *single* MWNTs and
       MWNT bundles.

    .. versionadded:: 0.2.8

    Parameters
    ----------
    Ch_list : :class:`python:list`, optional
        (:attr:`~SWNT.n`, :attr:`~SWNT.m`) for each `SWNT` wall in `MWNT`.
    Nwalls : int, optional
        Number of `SWNT` walls in `MWNT`.
    Lz : float, optional
        `MWNT` length in **Angstroms**.
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
    nx, ny : int, optional
        Number of repeat unit cells in the :math:`x, y` dimensions.
    vdw_radius : float, optional
        van der Waals radius of nanotube atoms
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
    bundle_packing : {'hcp', 'ccp'}, optional
        Packing arrangement of MWNT bundles.  If `bundle_packing` is `None`,
        then it will be determined by the `bundle_geometry` parameter if
        `bundle_geometry` is not `None`.  If both `bundle_packing` and
        `bundle_geometry` are `None`, then `bundle_packing` defaults to `hcp`.
    bundle_geometry : {'triangle', 'hexagon', 'square', 'rectangle'}, optional
    autogen : bool, optional
        if `True`, automatically call
        :meth:`~MWNTGenerator.generate`.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.generators import MWNTGenerator
    >>> MWNTGenerator(Nwalls=5, min_wall_diameter=10, Lz=50).save()

    The above command generated a 5 wall,  50 **Angstrom** long `MWNT`.
    The only constraints on the `MWNT` wall chiralities were the diameter
    constraints imposed by the `min_wall_diameter` parameter,
    which set the minimum wall diameter to 10 Angstroms,
    as well as the minimum wall-to-wall separation, which
    defaults to the van der Waals distance of 3.4 Angstroms.

    This `MWNT` chirality may be written as:

    :math:`\\mathbf{C}_{\\mathrm{h}} = (8,7)@(17,8)@(9,24)@(27,18)@(22,32)`

    Here's a colorful rendering of the generated `MWNT` structure:

    .. image:: /images/5wall_mwnt_(8,7)@(17,8)@(9,24)@(27,18)@(22,32)-04.png

    In the next example, we generate a `MWNT` bundle.

    >>> mwntbundle = MWNTGenerator(Nwalls=5, min_wall_diameter=5, Lz=50,
    ...                            bundle_geometry='hexagon')
    >>> mwntbundle.save()

    .. image:: /images/5wall_mwnt_(1,6)@(11,7)@(18,9)@(31,2)@(41,0)_hcp_7tube_hexagon-perspective_view-01.png

    """
    @classmethod
    def generate_fname(cls, Ch_list=None, Nwalls=None, Ntubes=None,
                       nx=None, ny=None, bundle_geometry=None,
                       bundle_packing=None, **kwargs):
        Nwalls = '{}wall_mwnt'.format(Nwalls)
        chiralities = '@'.join([str(Ch).replace(' ', '') for Ch in Ch_list])
        if nx == ny == 1 and bundle_geometry is None:
            fname = '_'.join((Nwalls, chiralities))
            return fname

        packing = '{}cp'.format(bundle_packing[0])
        Ntubes = '{}tube'.format(Ntubes)

        fname_wordlist = None
        if bundle_geometry is None:
            nx = ''.join(('{}'.format(nx), pluralize('cell', nx)))
            ny = ''.join(('{}'.format(ny), pluralize('cell', ny)))
            cells = 'x'.join((nx, ny))
            fname_wordlist = (Nwalls, chiralities, packing, cells)
        else:
            fname_wordlist = \
                (Nwalls, chiralities, packing, Ntubes, bundle_geometry)

        fname = '_'.join(fname_wordlist)
        return fname

    def save(self, fname=None, outpath=None, structure_format=None,
             center_centroid=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            fname = self.generate_fname(Ch_list=self.Ch_list,
                                        Nwalls=self.Nwalls,
                                        Ntubes=self.Ntubes,
                                        nx=self.nx, ny=self.ny,
                                        bundle_geometry=self.bundle_geometry,
                                        bundle_packing=self.bundle_packing)

        super().save(fname=fname, outpath=outpath,
                     structure_format=structure_format,
                     center_centroid=center_centroid, **kwargs)
