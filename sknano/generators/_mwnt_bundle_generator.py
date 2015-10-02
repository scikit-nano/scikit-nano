# -*- coding: utf-8 -*-
"""
==============================================================================
MWNT Bundle Generator (:mod:`sknano.generators._mwnt_bundle_generator`)
==============================================================================

.. currentmodule:: sknano.generators._mwnt_bundle_generator

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.core import pluralize
from sknano.structures import MWNTBundle
from ._mwnt_generator import MWNTGenerator
from ._nanotube_bundle_generator import NanotubeBundleGeneratorBase

__all__ = ['MWNTBundleGenerator']


class MWNTBundleGenerator(NanotubeBundleGeneratorBase, MWNTBundle,
                          MWNTGenerator):
    """Class for generating multi-walled nanotube bundles.

    .. versionadded:: 0.2.20

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
    nx, ny : int, optional
        Number of repeat unit cells in the :math:`x, y` dimensions.
    vdw_radius : float, optional
        van der Waals radius of nanotube atoms
    bundle_packing : {'hcp', 'ccp'}, optional
        Packing arrangement of MWNT bundles.  If `bundle_packing` is `None`,
        then it will be determined by the `bundle_geometry` parameter if
        `bundle_geometry` is not `None`.  If both `bundle_packing` and
        `bundle_geometry` are `None`, then `bundle_packing` defaults to `hcp`.
    bundle_geometry : {'triangle', 'hexagon', 'square', 'rectangle'}, optional
    autogen : bool, optional
        if `True`, automatically call
        :meth:`~MWNTBundleGenerator.generate`.
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

    >>> from sknano.generators import MWNTBundleGenerator
    >>> mwntbundle = MWNTBundleGenerator(Nwalls=5, min_wall_diameter=5, Lz=5,
    ...                                  bundle_geometry='hexagon')
    >>> mwntbundle.save()

    .. image:: /images/5wall_mwnt_(1,6)@(11,7)@(18,9)@(31,2)@(41,0)_hcp_7tube_hexagon-perspective_view-01.png

    """
    @classmethod
    def generate_fname(cls, Ch_list=None, Nwalls=None, Ntubes=None,
                       nx=None, ny=None, bundle_geometry=None,
                       bundle_packing=None, **kwargs):
        Nwalls = '{}wall_mwnt'.format(Nwalls)
        chiralities = '@'.join([str(Ch).replace(' ', '')
                                for Ch in Ch_list])
        packing = '{}cp'.format(bundle_packing[0])
        Ntubes = '{}tube'.format(Ntubes)

        fname_wordlist = None
        if bundle_geometry is None:
            nx = ''.join(('{}'.format(nx), pluralize('cell', nx)))
            ny = ''.join(('{}'.format(ny), pluralize('cell', ny)))
            cells = 'x'.join((nx, ny))

            if nx == 1 and ny == 1:
                fname_wordlist = (Nwalls, chiralities, cells)
            else:
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
