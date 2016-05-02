# -*- coding: utf-8 -*-
"""
===============================================================================
SWNT structure generator (:mod:`sknano.generators.swnt_generator`)
===============================================================================

.. currentmodule:: sknano.generators.swnt_generator

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
from sknano.core.structures import SWNT
from sknano.core.geometric_regions import Cuboid
from .base import NanoStructureGenerator, GeneratorMixin

from .nanotube_bundle_generator import NanotubeBundleGeneratorBase


__all__ = ['SWNTGeneratorBase', 'SWNTGenerator']


class SWNTGeneratorBase(GeneratorMixin, NanoStructureGenerator):
    """Base :class:`SWNT` generator class."""
    pass


class SWNTGenerator(NanotubeBundleGeneratorBase, SWNTGeneratorBase, SWNT):
    """Class for generating :class:`SWNT` structures.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of integers (i.e., *Ch = [(n, m)]) or
        2 integers (i.e., *Ch = [n, m]) specifying the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    n1, n2, n3 : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions.
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
    L : float, optional
        Length of nanotube in units of **Angstroms**.
        Overrides the `n3` value.

        .. versionadded:: 0.2.5

        .. versionchanged:: 0.4.0

           Changed units from nanometers to **Angstroms**

    tube_length : float, optional
        Length of nanotube in units of **Angstroms**.
        Overrides the `n3` value.

        .. deprecated:: 0.2.5
           Use `L` instead

    fix_L : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    vdw_radius : float, optional
        van der Waals radius of nanotube atoms

        .. versionadded:: 0.2.5

    bundle_packing : {'hcp', 'ccp'}, optional
        Packing arrangement of nanotubes bundles.
        If `bundle_packing` is `None`, then it will be determined by the
        `bundle_geometry` parameter if `bundle_geometry` is not `None`.
        If both `bundle_packing` and `bundle_geometry` are `None`, then
        `bundle_packing` defaults to `hcp`.

        .. versionadded:: 0.2.5

    bundle_geometry : {'triangle', 'hexagon', 'square', 'rectangle'}, optional
        Force a specific geometry on the nanotube bundle boundaries.

        .. versionadded:: 0.2.5

    autogen : bool, optional
        if `True`, automatically generate structure data.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------
    Using the `SWNTGenerator` class, you can generate structure
    data for single-walled nanotubes (SWNTs) or nanotube *bundles* with
    either cubic close packed (ccp) or hexagonal close packed (hcp)
    arrangement of nanotubes. The bundle packing arrangement is set using the
    `bundle_packing` parameter.

    You can also enforce a specific
    `bundle geometry` which will try and build the nanotube bundle such
    that it "fits inside" the boundaries of a specified geometric shape.
    This allows you to generate **hcp** bundles that are trianglar,
    hexagonal, or rectangular in shape, as some of the examples below
    illustrate.

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

    Next, let's generate an hcp bundle of
    :math:`C_{\\mathrm{h}} = (10, 5)` SWCNTs and cell count
    :math:`n_x=10, n_y=3, n_z=5`.

    >>> SWCNTbundle = SWNTGenerator(n=10, m=5, n1=10, n2=3, n3=5)
    >>> SWCNTbundle.save()

    The rendered structure looks like:

    .. image:: /images/1005_hcp_10cellsx3cellsx5cells-01.png

    Now let's generate a nice hexagon bundle, 3 tubes wide, with
    :math:`C_{\\mathrm{h}} = (6, 5)`.

    >>> SWCNTbundle = SWNTGenerator(n=6, m=5, n3=5, bundle_geometry='hexagon')
    >>> SWCNTbundle.save()

    which looks like:

    .. image:: /images/0605_hcp_7tube_hexagon-01.png

    Remember, all of the :meth:`~SWNTGenerator.save`
    methods allow you to rotate the structure data before saving:

    >>> SWCNTbundle.save(fname='0605_hcp_7tube_hexagon_rot-30deg.xyz',
    ...                  rot_axis='z', rotation_angle=30)

    .. image:: /images/0605_hcp_7tube_hexagon_rot-30deg-01.png

    Now, just because we can, let's make a big ass hexagon bundle with
    :math:`C_{\\mathrm{h}} = (10, 0)`.

    >>> big_hexagonal_bundle = SWNTGenerator(10, 0, n1=25, n2=25, n3=1,
    ...                                      bundle_geometry='hexagon')
    >>> big_hexagonal_bundle.save()

    Take a look at the 469 :math:`(10, 0)` unit cells in this big ass bundle!

    .. image:: /images/1000_hcp_469tube_hexagon-01.png

    Lastly, here's a look at a bundle generated with cubic close packed
    bundle arrangement:

    >>> SWCNTbundle = SWNTGenerator(10, 10, n1=3, n2=3, n3=5,
    ...                             bundle_packing='ccp')
    >>> SWCNTbundle.save()

    The rendered `ccp` structure looks like:

    .. image:: /images/1010_ccp_3cellsx3cellsx5cells-01.png

    """
    def finalize(self):
        """Finalize structure data by clipping region bounds if \
                :attr:`~SWNT.fix_L` is `True`."""
        if self.L0 is not None and self.fix_L:
            L_cutoff = self.L0 + self.bond
            region_bounds = Cuboid(pmin=Point([-np.inf, -np.inf, 0]),
                                   pmax=Point([np.inf, np.inf, L_cutoff]))
            self.atoms.clip_bounds(region_bounds)
        super().finalize()

    @classmethod
    def generate_fname(cls, n=None, m=None, n1=None, n2=None, n3=None,
                       fix_L=False, Ntubes=None, bundle_geometry=None,
                       bundle_packing=None, **kwargs):
        """Generate filename string."""
        fname = '{}{}'.format('{}'.format(n).zfill(2), '{}'.format(m).zfill(2))

        if n1 == n2 == 1 and bundle_geometry is None:
            n3_fmtstr = '{:.2f}' if fix_L else '{:.0f}'
            n3 = ''.join((n3_fmtstr.format(n3), pluralize('cell', n3)))
            return '_'.join((fname, n3))

        packing = '{}cp'.format(bundle_packing[0])
        Ntubes = '{}tube'.format(Ntubes)
        n1 = ''.join(('{}'.format(n1), pluralize('cell', n1)))
        n2 = ''.join(('{}'.format(n2), pluralize('cell', n2)))
        n3_fmtstr = '{:.2f}' if fix_L else '{:.0f}'
        n3 = ''.join((n3_fmtstr.format(n3), pluralize('cell', n3)))
        cells = 'x'.join((n1, n2, n3))

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
            fname = self.generate_fname(n=self.n, m=self.m,
                                        n1=self.n1, n2=self.n2, n3=self.n3,
                                        fix_L=self.fix_L,
                                        Ntubes=self.Ntubes,
                                        bundle_geometry=self.bundle_geometry,
                                        bundle_packing=self.bundle_packing)

        super().save(fname=fname, outpath=outpath,
                     structure_format=structure_format,
                     center_centroid=center_centroid, **kwargs)
