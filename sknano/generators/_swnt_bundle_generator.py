# -*- coding: utf-8 -*-
"""
==============================================================================
SWNT Bundle Generator (:mod:`sknano.generators._swnt_bundle_generator`)
==============================================================================

.. currentmodule:: sknano.generators._swnt_bundle_generator

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.core import pluralize
from sknano.structures import SWNTBundle
from ._swnt_generator import SWNTGenerator
from ._nanotube_bundle_generator import NanotubeBundleGeneratorBase

__all__ = ['SWNTBundleGenerator']


class SWNTBundleGenerator(NanotubeBundleGeneratorBase, SWNTGenerator,
                          SWNTBundle):
    """Class for generating nanotube bundles.

    .. versionadded:: 0.2.4

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of ints or 2 integers giving the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nx, ny, nz : int, optional
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

    Lz : float, optional
        length of bundle in :math:`z` dimension in **nanometers**.
        Overrides `nz` value.

        .. versionadded:: 0.2.5

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

    Using the `SWNTBundleGenerator` class, you can generate structure
    data for nanotube *bundles* with either cubic close packed (ccp) or
    hexagonal close packed (hcp) arrangement of nanotubes. The bundle
    packing arrangement is set using the `bundle_packing` parameter.

    You can also enforce a specific
    `bundle geometry` which will try and build the nanotube bundle such
    that it "fits inside" the boundaries of a specified geometric shape.
    This allows you to generate **hcp** bundles that are trianglar,
    hexagonal, or rectangular in shape, as some of the examples below
    illustrate.

    To start, let's generate an hcp bundle of
    :math:`C_{\\mathrm{h}} = (10, 5)` SWCNTs and cell count
    :math:`n_x=10, n_y=3, n_z=5`.

    >>> from sknano.generators import SWNTBundleGenerator
    >>> SWCNTbundle = SWNTBundleGenerator(n=10, m=5, nx=10, ny=3, nz=5)
    >>> SWCNTbundle.save()

    The rendered structure looks like:

    .. image:: /images/1005_hcp_10cellsx3cellsx5cells-01.png

    Now let's generate a nice hexagon bundle, 3 tubes wide, with
    :math:`C_{\\mathrm{h}} = (6, 5)`.

    >>> SWCNTbundle = SWNTBundleGenerator(n=6, m=5, nz=5,
    ...                                   bundle_geometry='hexagon')
    >>> SWCNTbundle.save()

    which looks like:

    .. image:: /images/0605_hcp_7tube_hexagon-01.png

    Remember, all of the :meth:`~SWNTBundleGenerator.save`
    methods allow you to rotate the structure data before saving:

    >>> SWCNTbundle.save(fname='0605_hcp_7tube_hexagon_rot-30deg.xyz',
    ...                  rot_axis='z', rotation_angle=30)

    .. image:: /images/0605_hcp_7tube_hexagon_rot-30deg-01.png

    Now, just because we can, let's make a big ass hexagon bundle with
    :math:`C_{\\mathrm{h}} = (10, 0)`.

    >>> big_hexagonal_bundle = SWNTBundleGenerator(10, 0, nx=25, ny=25, nz=1,
    ...                                            bundle_geometry='hexagon')
    >>> big_hexagonal_bundle.save()

    Take a look at the 469 :math:`(10, 0)` unit cells in this big ass bundle!

    .. image:: /images/1000_hcp_469tube_hexagon-01.png

    Lastly, here's a look at a bundle generated with cubic close packed
    bundle arrangement:

    >>> SWCNTbundle = SWNTBundleGenerator(10, 10, nx=3, ny=3, nz=5,
    ...                                   bundle_packing='ccp')
    >>> SWCNTbundle.save()

    The rendered `ccp` structure looks like:

    .. image:: /images/1010_ccp_3cellsx3cellsx5cells-01.png

    """
    @classmethod
    def generate_fname(cls, n=None, m=None, nx=None, ny=None, nz=None,
                       fix_Lz=False, Ntubes=None, bundle_geometry=None,
                       bundle_packing=None, **kwargs):
        chirality = '{}{}'.format('{}'.format(n).zfill(2),
                                  '{}'.format(m).zfill(2))
        packing = '{}cp'.format(bundle_packing[0])
        Ntube = '{}tube'.format(Ntubes)

        fname_wordlist = None
        if bundle_geometry is None:
            nx = ''.join(('{}'.format(nx),
                         pluralize('cell', nx)))
            ny = ''.join(('{}'.format(ny),
                         pluralize('cell', ny)))
            nz_fmtstr = '{:.2f}' if fix_Lz else '{:.0f}'
            nz = ''.join((nz_fmtstr.format(nz), pluralize('cell', nz)))
            cells = 'x'.join((nx, ny, nz))
            fname_wordlist = (chirality, packing, cells)
        else:
            fname_wordlist = \
                (chirality, packing, Ntube, bundle_geometry)

        fname = '_'.join(fname_wordlist)
        return fname

    def save(self, fname=None, outpath=None, structure_format=None,
             center_CM=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            fname = self.generate_fname(n=self.n, m=self.m,
                                        nx=self.nx, ny=self.ny, nz=self.nz,
                                        fix_Lz=self.fix_Lz,
                                        Ntubes=self.Ntubes,
                                        bundle_geometry=self.bundle_geometry,
                                        bundle_packing=self.bundle_packing)

        super().save(fname=fname, outpath=outpath,
                     structure_format=structure_format,
                     center_CM=center_CM, **kwargs)
