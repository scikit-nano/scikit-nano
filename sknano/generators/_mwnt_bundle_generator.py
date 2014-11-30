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
from ._mixins import NanotubeBundleGeneratorMixin

__all__ = ['MWNTBundleGenerator']


class MWNTBundleGenerator(NanotubeBundleGeneratorMixin, MWNTBundle,
                          MWNTGenerator):
    """Class for generating multi-walled nanotube bundles.

    .. versionadded:: 0.2.20

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
    Lx, Ly, Lz : float, optional
        length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.
    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.
    add_outer_shells : bool, optional
        Build the MWNT by adding outer shells

        .. versionadded:: 0.2.23

    add_inner_shells : bool, optional
        Build the MWNT by adding inner shells

        .. versionadded:: 0.2.23

    max_shells : int, optional
        Maximum number of shells per MWNT.
    max_shell_diameter : float, optional
        Maximum shell diameter, in units of **Angstroms**.
    min_shells : int, optional
        Minimum number of shells per MWNT.
    min_shell_diameter : float, optional
        Minimum shell diameter, in units of **Angstroms**.
    new_shell_type : {None, 'AC', 'ZZ', 'achiral', 'chiral'}, optional
        If `None`, the chiralities of the new shells are constrained only
        by their diameter and will be chosen randomly if more than one
        candidate chirality exists. If not `None`, then the
        `new_shell_type` will be added as a constraint.
    shell_spacing : float, optional
        Shell spacing in units of **Angstroms**. Default
        value is the van der Waals interaction distance of 3.4 Angstroms.
    vdw_spacing : float, optional
        van der Waals distance between nearest neighbor tubes
    bundle_packing : {'hcp', 'ccp'}, optional
        Packing arrangement of MWNT bundles.  If `bundle_packing` is `None`,
        then it will be determined by the `bundle_geometry` parameter if
        `bundle_geometry` is not `None`.  If both `bundle_packing` and
        `bundle_geometry` are `None`, then `bundle_packing` defaults to `hcp`.
    bundle_geometry : {'triangle', 'hexagon', 'square', 'rectangle'}, optional
    autogen : bool, optional
        if `True`, automatically call
        :meth:`~MWNTGenerator.generate_unit_cell`,
        followed by :meth:`~MWNTGenerator.generate_structure_data`.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.generators import MWNTBundleGenerator
    >>> mwntbundle = MWNTBundleGenerator(max_shells=5, Lz=1.0, fix_Lz=True,
    ...                                  bundle_geometry='hexagon')
    >>> mwntbundle.save_data()

    .. image:: /images/5shell_mwnt_4040_outer_Ch_hcp_7tube_hexagon-01.png

    """
    def __init__(self, autogen=True, **kwargs):

        super(MWNTBundleGenerator, self).__init__(autogen=False, **kwargs)

        if autogen:
            super(MWNTBundleGenerator, self).generate_unit_cell()
            self.generate_structure_data()

    def generate_structure_data(self):
        """Generate structure data."""
        super(MWNTBundleGenerator, self).generate_structure_data()
        super(MWNTBundleGenerator, self).generate_bundle()

    def save_data(self, fname=None, outpath=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, anchor_point=None,
                  deg2rad=True, center_CM=True, savecopy=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save_data` method
        for documentation.

        """
        if fname is None:
            Nshells = '{}shell_mwnt'.format(self.Nshells)
            chiralities = '@'.join([str(Ch).replace(' ', '')
                                    for Ch in self.Ch])
            packing = '{}cp'.format(self.bundle_packing[0])
            Ntubes = '{}tube'.format(self.Ntubes)

            fname_wordlist = None
            if self.bundle_geometry is None:
                nx = ''.join(('{}'.format(self.nx),
                             pluralize('cell', self.nx)))
                ny = ''.join(('{}'.format(self.ny),
                             pluralize('cell', self.ny)))
                cells = 'x'.join((nx, ny))
                #if self._assert_integer_nz:
                #    nz = ''.join(('{}'.format(self.nz),
                #                  pluralize('cell', self.nz)))
                #else:
                #    nz = ''.join(('{:.2f}'.format(self.nz),
                #                  pluralize('cell', self.nz)))
                #cells = 'x'.join((nx, ny, nz))

                if self.nx == 1 and self.ny == 1:
                    fname_wordlist = (Nshells, chiralities, cells)
                else:
                    fname_wordlist = (Nshells, chiralities, packing, cells)
            else:
                fname_wordlist = \
                    (Nshells, chiralities, packing, Ntubes,
                     self.bundle_geometry)

            fname = '_'.join(fname_wordlist)

        super(MWNTBundleGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            anchor_point=anchor_point, deg2rad=deg2rad, center_CM=center_CM,
            savecopy=savecopy, **kwargs)
