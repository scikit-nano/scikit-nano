# -*- coding: utf-8 -*-
"""
===============================================================================
MWNT structure generator (:mod:`sknano.generators.mwnt_generator`)
===============================================================================

.. currentmodule:: sknano.generators.mwnt_generator

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

# import copy

import numpy as np

from sknano.core import pluralize
from sknano.core.math import Point
from sknano.core.structures import MWNT
from sknano.core.geometric_regions import Cuboid
from .base import GeneratorMixin, NanoStructureGenerator
from .nanotube_bundle_generator import NanotubeBundleGeneratorBase
# from .swnt_generator import SWNTGenerator

__all__ = ['MWNTGeneratorBase', 'MWNTGenerator']


class MWNTGeneratorBase(GeneratorMixin, NanoStructureGenerator):
    """Mixin :class:`~sknano.core.structures.MWNT` structure generator."""
    pass


class MWNTGenerator(NanotubeBundleGeneratorBase, MWNTGeneratorBase, MWNT):
    """Class for generating multi-walled nanotubes (MWNT).

    .. versionchanged:: 0.4.0

       `MWNTGenerator` now generates both *single* MWNTs and
       MWNT bundles.

    .. versionadded:: 0.2.8

    Parameters
    ----------
    Ch_list : :class:`python:list`, optional
        (:attr:`~SWNT.n`, :attr:`~SWNT.m`) for each `SWNT` wall in `MWNT`.
    Nwalls : int, optional
        Number of `SWNT` walls in `MWNT`.
    L : float, optional
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
    n1, n2 : int, optional
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
    >>> MWNTGenerator(Nwalls=5, min_wall_diameter=10, L=50).save()

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

    >>> mwntbundle = MWNTGenerator(Nwalls=5, min_wall_diameter=5, L=50,
    ...                            bundle_geometry='hexagon')
    >>> mwntbundle.save()

    .. image:: /images/5wall_mwnt_(1,6)@(11,7)@(18,9)@(31,2)@(41,0)_hcp_7tube_hexagon-perspective_view-01.png

    """
    def finalize(self):
        """Finalize structure data by clipping region bounds if \
                :attr:`~MWNT.fix_L` is `True`."""
        # print('finalizing structure data')
        # print('Natoms: {}'.format(self.atoms.Natoms))
        if self.fix_L:
            L_cutoff = self.L + self.bond
            # print('clipping_bounds')
            # print(self.atoms.coordinates_bounding_box)
            # print('L_cutoff: {}'.format(L_cutoff))
            region_bounds = Cuboid(pmin=Point([-np.inf, -np.inf, 0]),
                                   pmax=Point([np.inf, np.inf, L_cutoff]))
            self.atoms.clip_bounds(region_bounds)
        # print('Natoms: {}'.format(self.atoms.Natoms))
        super().finalize()

    @classmethod
    def generate_fname(cls, Ch_list=None, Nwalls=None, Ntubes=None, L=None,
                       n1=None, n2=None, n3=None, bundle_geometry=None,
                       bundle_packing=None, **kwargs):
        """Generate a filename string."""
        # Nwalls = '{}wall_mwnt'.format(Nwalls)
        fname = '@'.join([str(Ch).replace(' ', '') for Ch in Ch_list])
        if n1 == n2 == 1 and bundle_geometry is None:
            return fname
            # fname = '_'.join((Nwalls, chiralities))
            # return fname

        packing = '{}cp'.format(bundle_packing[0])
        Ntubes = '{}tube'.format(Ntubes)

        n1 = ''.join(('{}'.format(n1), pluralize('cell', n1)))
        n2 = ''.join(('{}'.format(n2), pluralize('cell', n2)))
        cells = 'x'.join((n1, n2))

        if n3 is not None:
            n3 = ''.join(('{}'.format(n3), pluralize('cell', n3)))
            cells = 'x'.join((cells, n3))
        else:
            cells = 'x'.join((cells, '{:.1f}Å'.format(L)))

        fname = '_'.join((fname, cells, packing))

        if bundle_geometry is not None:
            fname = '_'.join((fname, Ntubes, bundle_geometry))
        return fname

    def save(self, fname=None, outpath=None, structure_format=None,
             center_centroid=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            fname = self.generate_fname(Ch_list=self.Ch_list,
                                        Nwalls=self.Nwalls, Ntubes=self.Ntubes,
                                        L=self.L,
                                        n1=self.n1, n2=self.n2, n3=self.n3,
                                        bundle_geometry=self.bundle_geometry,
                                        bundle_packing=self.bundle_packing)

        super().save(fname=fname, outpath=outpath,
                     structure_format=structure_format,
                     center_centroid=center_centroid, **kwargs)
