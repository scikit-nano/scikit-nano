# -*- coding: utf-8 -*-
"""
==============================================================================
Bundle Generators (:mod:`sknano.generators._nanotube_bundle_generators`)
==============================================================================

.. currentmodule:: sknano.generators._nanotube_bundle_generators

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import copy
#import itertools

import numpy as np

from ..core import plural_word_check
from ..core.atoms import Atoms
from ..structures import SWNTBundle, MWNTBundle
from ._nanotube_generators import SWNTGenerator, MWNTGenerator

__all__ = ['SWNTBundleGenerator', 'MWNTBundleGenerator']


class NanotubeBundleGeneratorMixin(object):

    def generate_hexagonal_bundle(self):
        nrows = max(self._nx, self._ny, 3)
        if nrows % 2 != 1:
            nrows += 1

        ntubes_per_end_rows = int((nrows + 1) / 2)

        row = 0
        ntubes_per_row = nrows
        while ntubes_per_row >= ntubes_per_end_rows:
            if row == 0:
                for n in xrange(ntubes_per_row):
                    atomsobj = Atoms(atoms=self._atomsobj0, deepcopy=True)
                    atomsobj.center_CM()
                    dr = n * self._r1
                    atomsobj.translate(dr)
                    self._structure_atoms.extend(atomsobj)
                    self._Ntubes += 1
            else:
                for nx in xrange(ntubes_per_row):
                    for ny in (-row, row):
                        atomsobj = Atoms(atoms=self._atomsobj0, deepcopy=True)
                        atomsobj.center_CM()
                        dr = np.zeros(3)
                        dr[0] = abs(ny * self._r2[0])
                        dr[1] = ny * self._r2[1]
                        dr = nx * self._r1 + dr
                        atomsobj.translate(dr)
                        self._structure_atoms.extend(atomsobj)
                        self._Ntubes += 1
            row += 1
            ntubes_per_row = nrows - row

    def generate_rectangular_bundle(self):
        Lx = 10 * self._Lx
        for nx in xrange(self._nx):
            for ny in xrange(self._ny):
                atomsobj = Atoms(atoms=self._atomsobj0, deepcopy=True)
                atomsobj.center_CM()
                dr = nx * self._r1 + ny * self._r2
                while dr[0] < 0:
                    dr[0] += Lx
                atomsobj.translate(dr)
                self._structure_atoms.extend(atomsobj)
                self._Ntubes += 1

    def generate_square_bundle(self):
        pass

    def generate_triangular_bundle(self):
        pass

    def generate_bundle(self):
        self._Ntubes = 0
        self._atomsobj0 = copy.deepcopy(self._structure_atoms)
        self._structure_atoms = Atoms()

        if self._bundle_geometry == 'hexagon':
            self.generate_hexagonal_bundle()
        elif self._bundle_geometry == 'rectangle':
            self.generate_rectangular_bundle()
        elif self._bundle_geometry == 'square':
            self.generate_square_bundle()
        elif self._bundle_geometry == 'triangle':
            self.generate_triangular_bundle()
        else:
            for nx in xrange(self._nx):
                for ny in xrange(self._ny):
                    atomsobj = Atoms(atoms=self._atomsobj0, deepcopy=True)
                    atomsobj.center_CM()
                    dr = nx * self._r1 + ny * self._r2
                    atomsobj.translate(dr)
                    self._structure_atoms.extend(atomsobj)
                    self._Ntubes += 1


class SWNTBundleGenerator(NanotubeBundleGeneratorMixin, SWNTBundle,
                          SWNTGenerator):
    u"""Class for generating nanotube bundles.

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

    Lx, Ly, Lz : float, optional
        length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.

        .. versionadded:: 0.2.5

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    autogen : bool, optional
        if `True`, automatically call
        :meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :meth:`~NanotubeGenerator.generate_structure_data`.
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
    ...                                       bundle_geometry='hexagon')
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
    ...                                         bundle_geometry='hexagon')
    >>> BIGASSHEXABUN.save_data()

    Take a look at the 469 :math:`(10, 0)` unit cells in this big ass bundle!

    .. image:: /images/1000_hcp_469tube_hexagon-01.png

    Lastly, here's a look at a bundle generated with cubic close packed
    bundle arrangement:

    >>> SWCNTbundle = SWNTBundleGenerator(10, 10, nx=3, ny=3, nz=5,
    ...                                       bundle_packing='ccp')
    >>> SWCNTbundle.save_data()

    The rendered `ccp` structure looks like:

    .. image:: /images/1010_ccp_3cellsx3cellsx5cells-01.png

    """

    def __init__(self, autogen=True, **kwargs):

        super(SWNTBundleGenerator, self).__init__(autogen=False, **kwargs)

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

        See :meth:`~sknano.generators.GeneratorMixin.save_data` method
        for documentation.

        """
        if fname is None:
            chirality = '{}{}'.format('{}'.format(self._n).zfill(2),
                                      '{}'.format(self._m).zfill(2))
            packing = '{}cp'.format(self._bundle_packing[0])
            #Ntubes = ''.join(('{}'.format(self._Ntubes),
            #                  plural_word_check('tube', self._Ntubes)))
            Ntube = '{}tube'.format(self._Ntubes)

            fname_wordlist = None
            if self._bundle_geometry is None:
                nx = ''.join(('{}'.format(self._nx),
                             plural_word_check('cell', self._nx)))
                ny = ''.join(('{}'.format(self._ny),
                             plural_word_check('cell', self._ny)))
                if self._assume_integer_unit_cells:
                    nz = ''.join(('{}'.format(self._nz),
                                  plural_word_check('cell', self._nz)))
                else:
                    nz = ''.join(('{:.2f}'.format(self._nz),
                                  plural_word_check('cell', self._nz)))
                cells = 'x'.join((nx, ny, nz))
                fname_wordlist = (chirality, packing, cells)
            else:
                fname_wordlist = \
                    (chirality, packing, Ntube, self._bundle_geometry)

            fname = '_'.join(fname_wordlist)

        super(SWNTBundleGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            anchor_point=anchor_point, deg2rad=deg2rad, center_CM=center_CM,
            savecopy=savecopy, **kwargs)


class MWNTBundleGenerator(NanotubeBundleGeneratorMixin, MWNTBundle,
                          MWNTGenerator):
    u"""Class for generating multi-walled nanotube bundles.

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
    >>> mwntbundle = MWNTBundleGenerator(n=40, m=40, max_shells=5,
    ...                                  Lz=1.0, fix_Lz=True,
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

        See :meth:`~sknano.generators.GeneratorMixin.save_data` method
        for documentation.

        """
        if fname is None:
            Nshells = '{}shell_mwnt'.format(self._Nshells_per_tube)
            chirality = '{}{}_{}_Ch'.format('{}'.format(self._n).zfill(2),
                                            '{}'.format(self._m).zfill(2),
                                            self._starting_shell_position)
            packing = '{}cp'.format(self._bundle_packing[0])
            Ntubes = '{}tube'.format(self._Ntubes)

            fname_wordlist = None
            if self._bundle_geometry is None:
                nx = ''.join(('{}'.format(self._nx),
                             plural_word_check('cell', self._nx)))
                ny = ''.join(('{}'.format(self._ny),
                             plural_word_check('cell', self._ny)))
                if self._assume_integer_unit_cells:
                    nz = ''.join(('{}'.format(self._nz),
                                  plural_word_check('cell', self._nz)))
                else:
                    nz = ''.join(('{:.2f}'.format(self._nz),
                                  plural_word_check('cell', self._nz)))
                cells = 'x'.join((nx, ny, nz))

                if self._nx == 1 and self._ny == 1:
                    fname_wordlist = (Nshells, chirality, cells)
                else:
                    fname_wordlist = (Nshells, chirality, packing, cells)
            else:
                fname_wordlist = \
                    (Nshells, chirality, packing, Ntubes,
                     self._bundle_geometry)

            fname = '_'.join(fname_wordlist)

        super(MWNTBundleGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            anchor_point=anchor_point, deg2rad=deg2rad, center_CM=center_CM,
            savecopy=savecopy, **kwargs)
