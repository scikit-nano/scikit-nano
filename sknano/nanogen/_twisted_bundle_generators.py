# -*- coding: utf-8 -*-
"""
==============================================================================
Twisted nanotube bundles (:mod:`sknano.nanogen._twisted_bundle_generators`)
==============================================================================

.. currentmodule:: sknano.nanogen._twisted_bundle_generators

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#import copy
#import itertools
#import sys

#import numpy as np

from ..tools.refdata import CCbond

from ._nanotube_bundle_generators import NanotubeBundleGenerator

__all__ = ['TwistedBundleGenerator']


class TwistedBundleGenerator(NanotubeBundleGenerator):
    u"""Class for generating twisted nanotube bundles.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :py:class:`~sknano.chemistry.Atoms` 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    vdw_spacing : float, optional
        van der Waals distance between nearest neighbor tubes
    bundle_packing : {None, 'hcp', 'hexagonal', 'ccp', 'cubic'}, optional
        Packing arrangement of nanotubes bundles.
        If `bundle_packing` is `None`, then it will be determined by the
        `bundle_geometry` parameter if `bundle_geometry` is not `None`.
        If both `bundle_packing` and `bundle_geometry` are `None`, then
        `bundle_packing` defaults to `hexagonal`.
    bundle_geometry : {None, 'triangle', 'hexagon', 'square', 'rectangle',
                       'rhombus', 'rhomboid'}, optional
        Force a specific geometry on the nanotube bundle boundaries.
    Lx, Ly, Lz : float, optional
        length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.
    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.
    autogen : bool, optional
        if `True`, automatically call
        :py:meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by
        :py:meth:`~NanotubeBundleGenerator.generate_structure_data`.
    verbose : bool, optional
        if `True`, show verbose output

    """

    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond, vdw_spacing=3.4,
                 bundle_packing=None, bundle_geometry=None, Lx=None, Ly=None,
                 Lz=None, fix_Lz=False, autogen=True, verbose=False):

        raise RuntimeError('This class is not yet implemented.')

        super(TwistedBundleGenerator, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz, element1=element1,
            element2=element2, Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz, bond=bond,
            autogen=False, verbose=verbose)

        if autogen:
            super(TwistedBundleGenerator, self).generate_unit_cell()
            super(TwistedBundleGenerator, self).generate_structure_data()
            self.twist_structure_data()

    def twist_structure_data(self):
        pass
