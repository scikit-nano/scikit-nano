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
from ._mixins import NanotubeBundleGeneratorMixin

__all__ = ['SWNTBundleGenerator']


class SWNTBundleGenerator(NanotubeBundleGeneratorMixin, SWNTBundle,
                          SWNTGenerator):
    """Class for generating nanotube bundles.

    .. versionadded:: 0.2.4

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.Atom` 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    vdw_spacing : float, optional
        van der Waals distance between nearest neighbor tubes

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
    >>> SWCNTbundle.save_data()

    The rendered structure looks like:

    .. image:: /images/1005_hcp_10cellsx3cellsx5cells-01.png

    Now let's generate a nice hexagon bundle, 3 tubes wide, with
    :math:`C_{\\mathrm{h}} = (6, 5)`.

    >>> SWCNTbundle = SWNTBundleGenerator(n=6, m=5, nz=5,
    ...                                   bundle_geometry='hexagon')
    >>> SWCNTbundle.save_data()

    which looks like:

    .. image:: /images/0605_hcp_7tube_hexagon-01.png

    Remember, all of the :meth:`~SWNTBundleGenerator.save_data`
    methods allow you to rotate the structure data before saving:

    >>> SWCNTbundle.save_data(fname='0605_hcp_7tube_hexagon_rot-30deg.xyz',
    ...                       rot_axis='z', rotation_angle=30)

    .. image:: /images/0605_hcp_7tube_hexagon_rot-30deg-01.png

    Now, just because we can, let's make a big ass hexagon bundle with
    :math:`C_{\\mathrm{h}} = (10, 0)`.

    >>> BIGASSHEXABUN = SWNTBundleGenerator(10, 0, nx=25, ny=25, nz=1,
    ...                                     bundle_geometry='hexagon')
    >>> BIGASSHEXABUN.save_data()

    Take a look at the 469 :math:`(10, 0)` unit cells in this big ass bundle!

    .. image:: /images/1000_hcp_469tube_hexagon-01.png

    Lastly, here's a look at a bundle generated with cubic close packed
    bundle arrangement:

    >>> SWCNTbundle = SWNTBundleGenerator(10, 10, nx=3, ny=3, nz=5,
    ...                                   bundle_packing='ccp')
    >>> SWCNTbundle.save_data()

    The rendered `ccp` structure looks like:

    .. image:: /images/1010_ccp_3cellsx3cellsx5cells-01.png

    """

    def __init__(self, *Ch, autogen=True, **kwargs):

        super(SWNTBundleGenerator, self).__init__(*Ch, autogen=False, **kwargs)

        if autogen:
            super(SWNTBundleGenerator, self).generate_unit_cell()
            self.generate_structure_data()

    def generate_structure_data(self):
        """Generate structure data."""
        super(SWNTBundleGenerator, self).generate_structure_data()
        super(SWNTBundleGenerator, self).generate_bundle()

    def save_data(self, fname=None, outpath=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, anchor_point=None,
                  deg2rad=True, center_CM=True, savecopy=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save_data` method
        for documentation.

        """
        if fname is None:
            chirality = '{}{}'.format('{}'.format(self.n).zfill(2),
                                      '{}'.format(self.m).zfill(2))
            packing = '{}cp'.format(self.bundle_packing[0])
            #Ntubes = ''.join(('{}'.format(self._Ntubes),
            #                  pluralize('tube', self._Ntubes)))
            Ntube = '{}tube'.format(self.Ntubes)

            fname_wordlist = None
            if self.bundle_geometry is None:
                nx = ''.join(('{}'.format(self.nx),
                             pluralize('cell', self.nx)))
                ny = ''.join(('{}'.format(self.ny),
                             pluralize('cell', self.ny)))
                if self._assert_integer_nz:
                    nz = ''.join(('{}'.format(self.nz),
                                  pluralize('cell', self.nz)))
                else:
                    nz = ''.join(('{:.2f}'.format(self.nz),
                                  pluralize('cell', self.nz)))
                cells = 'x'.join((nx, ny, nz))
                fname_wordlist = (chirality, packing, cells)
            else:
                fname_wordlist = \
                    (chirality, packing, Ntube, self.bundle_geometry)

            fname = '_'.join(fname_wordlist)

        super(SWNTBundleGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            anchor_point=anchor_point, deg2rad=deg2rad, center_CM=center_CM,
            savecopy=savecopy, **kwargs)
