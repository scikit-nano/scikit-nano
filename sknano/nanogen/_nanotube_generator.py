# -*- coding: utf-8 -*-
"""
=====================================================================
Nanotube structure tools (:mod:`sknano.nanogen._nanotube_generator`)
=====================================================================

.. currentmodule:: sknano.nanogen._nanotube_generator

.. todo::

   Add methods to perform fractional translation and cartesian translation
   before structure generation.

.. todo::

   Handle different units and perform unit conversions for output coordinates.

.. todo::

   Provide options for setting gutter (van der Waals separation) value

.. todo::

   Consider replacing coordinate arrays with class attributes to make
   code and math operations more readable

"""
from __future__ import division, print_function, absolute_import
__docformat__ = 'restructuredtext'

import copy
#import itertools
#import sys
import warnings
warnings.filterwarnings('ignore')  # to suppress the Pint UnicodeWarning

#from pint import UnitRegistry
#ureg = UnitRegistry()
#Qty = ureg.Quantity

import numpy as np

from pkshared.tools.arrayfuncs import rotation_matrix
from pkshared.tools.strfuncs import plural_word_check
from pkshared.tools.refdata import CCbond

from ..chemistry import Atom, Atoms
from ..structure_io import DATAWriter, XYZWriter, default_structure_format, \
    supported_structure_formats

from ._nanotube import Nanotube, NanotubeBundle

__all__ = ['NanotubeGeneratorError',
           'NanotubeGenerator',
           'NanotubeBundleGenerator',
           'MWNTGenerator']


class NanotubeGeneratorError(Exception):
    """Base class for NanotubeGenerator exceptions."""
    pass


class NanotubeGenerator(Nanotube):
    u"""Class for generating nanotube structures.

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nz : int, optional
        Number of repeat unit cells in the :math:`z` direction, along
        the *length* of the nanotube.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    Lz : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the ``nz`` value.

        .. versionadded:: 0.2.5

    tube_length : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the ``nz`` value.

        .. deprecated:: 0.2.5
           Use ``Lz`` instead

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If ``True``, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    autogen : bool, optional
        if ``True``, automatically call
        :py:meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :py:meth:`~NanotubeGenerator.generate_structure_data`.
    verbose : bool, optional
        if ``True``, show verbose output

    Examples
    --------
    First, load the :py:class:`~sknano.nanogen.NanoGenerator` class.

    >>> from sknano.nanogen import NanotubeGenerator

    Now let's generate a :math:`\\mathbf{C}_{\\mathrm{h}} = (10, 5)`
    SWCNT unit cell.

    >>> nt = NanotubeGenerator(n=10, m=5)
    >>> nt.save_data(fname='10,5_unit_cell.xyz')

    The rendered structure looks like (orhographic view):

    .. image:: /images/10,5_unit_cell_orthographic_view.png

    and the perspective view:

    .. image:: /images/10,5_unit_cell_perspective_view.png

    """

    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C',
                 bond=CCbond, Lx=None, Ly=None, Lz=None,
                 tube_length=None, fix_Lz=False,
                 autogen=True, verbose=False):

        if tube_length is not None and Lz is None:
            Lz = tube_length

        super(NanotubeGenerator, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz,
            element1=element1, element2=element2,
            bond=bond, Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz,
            with_units=False, verbose=verbose)

        self._fname = None
        self.unit_cell = None
        self.structure_atoms = None

        if autogen:
            self.generate_unit_cell()
            self.generate_structure_data()

    @property
    def fname(self):
        """Structure file name."""
        return self._fname

    def generate_unit_cell(self):
        """Generate the nanotube unit cell."""
        eps = 0.01
        n = self._n
        m = self._m
        bond = self._bond
        M = self._M
        T = self._T
        N = self._N
        rt = self._rt
        e1 = self._element1
        e2 = self._element2
        verbose = self._verbose

        aCh = Nanotube.compute_chiral_angle(n=n, m=m, rad2deg=False)

        tau = M * T / N
        dtau = bond * np.sin(np.pi / 6 - aCh)

        psi = 2 * np.pi / N
        dpsi = bond * np.cos(np.pi / 6 - aCh) / rt

        if verbose:
            print('dpsi: {}'.format(dpsi))
            print('dtau: {}\n'.format(dtau))

        self.unit_cell = Atoms()

        for i in xrange(1, N + 1):
            x1 = rt * np.cos(i * psi)
            y1 = rt * np.sin(i * psi)
            z1 = i * tau

            while z1 > T - eps:
                z1 -= T

            atom1 = Atom(e1, x=x1, y=y1, z=z1)
            atom1.rezero_coords()

            if verbose:
                print('Basis Atom 1:\n{}'.format(atom1))

            self.unit_cell.append(atom1)

            x2 = rt * np.cos(i * psi + dpsi)
            y2 = rt * np.sin(i * psi + dpsi)
            z2 = i * tau - dtau
            while z2 > T - eps:
                z2 -= T

            atom2 = Atom(e2, x=x2, y=y2, z=z2)
            atom2.rezero_coords()

            if verbose:
                print('Basis Atom 2:\n{}'.format(atom2))

            self.unit_cell.append(atom2)

    def generate_structure_data(self):
        """Generate structure data."""
        self.structure_atoms = []
        for nz in xrange(int(np.ceil(self._nz))):
            dr = np.array([0.0, 0.0, nz * self.T])
            for uc_atom in self.unit_cell.atoms:
                nt_atom = Atom(uc_atom.symbol)
                nt_atom.r = uc_atom.r + dr
                self.structure_atoms.append(nt_atom)

    def save_data(self, fname=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True):
        """Save structure data.

        Parameters
        ----------
        fname : {None, str}, optional
            file name string
        structure_format : {None, str}, optional
            chemical file format of saved structure data. Must be one of:

                - xyz
                - data

            If ``None``, then guess based on ``fname`` file extension or
            default to ``xyz`` format.
        rotation_angle : {None, float}, optional
            Angle of rotation
        rot_axis : {'x', 'y', 'z'}, optional
            Rotation axis
        deg2rad : bool, optional
            Convert ``rotation_angle`` from degrees to radians.
        center_CM : bool, optional
            Center center-of-mass on origin.

        """
        if (fname is None and structure_format not in
                supported_structure_formats) or \
                (fname is not None and not
                    fname.endswith(supported_structure_formats) and
                    structure_format not in supported_structure_formats):
            structure_format = default_structure_format

        if fname is None:
            chirality = '{}{}r'.format('{}'.format(self._n).zfill(2),
                                       '{}'.format(self._m).zfill(2))
            if self._assume_integer_unit_cells:
                nz = ''.join(('{}'.format(self._nz),
                              plural_word_check('cell', self._nz)))
            else:
                nz = ''.join(('{:.2f}'.format(self._nz),
                              plural_word_check('cell', self._nz)))
            fname_wordlist = (chirality, nz)
            fname = '_'.join(fname_wordlist)
            fname += '.' + structure_format
        else:
            if fname.endswith(supported_structure_formats) and \
                    structure_format is None:
                for ext in supported_structure_formats:
                    if fname.endswith(ext):
                        structure_format = ext
                        break
            else:
                if structure_format is None or \
                        structure_format not in supported_structure_formats:
                    structure_format = default_structure_format

        #structure_atoms = list(itertools.chain(*self.structure_atoms))
        structure_atoms = None
        if isinstance(self.structure_atoms, list):
            structure_atoms = Atoms(self.structure_atoms)
        elif isinstance(self.structure_atoms, Atoms):
            structure_atoms = self.structure_atoms

        if center_CM:
            structure_atoms.center_CM()

        if self._L0 is not None and self._fix_Lz:
            structure_atoms.clip_bounds(abs_limit=(10 * self._L0 + 0.5) / 2,
                                        r_indices=[2])

        if rotation_angle is not None:
            R_matrix = rotation_matrix(rotation_angle,
                                       rot_axis=rot_axis,
                                       deg2rad=deg2rad)
            structure_atoms.rotate(R_matrix)

        if structure_format == 'data':
            DATAWriter.write(fname=fname, atoms=structure_atoms)
        else:
            XYZWriter.write(fname=fname, atoms=structure_atoms)

        self._fname = fname


class NanotubeBundleGenerator(NanotubeGenerator, NanotubeBundle):
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
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    vdw_spacing : float, optional
        van der Waals distance between nearest neighbor tubes

        .. versionadded:: 0.2.5

    bundle_packing : {None, 'hexagonal', 'cubic'}, optional
        close packing arrangement of bundles

        .. versionadded:: 0.2.5

    bundle_geometry : {None, 'triangle', 'hexagon', 'square', 'rectangle',
                       'rhombus', 'rhomboid'}, optional

        .. versionadded:: 0.2.5

    Lx, Ly, Lz : float, optional
        length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.

        .. versionadded:: 0.2.5

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If ``True``, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    autogen : bool, optional
        if ``True``, automatically call
        :py:meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :py:meth:`~NanotubeGenerator.generate_structure_data`.
    verbose : bool, optional
        if ``True``, show verbose output

    Examples
    --------

    Using the :py:class:`NanotubeBundleGenerator` class, you can
    generate cubic close packed (ccp) or hexagonal close packed
    bundles arrangements. In general, specifying **ccp** bundling will
    generate rectangular bundles (square bundles if :math:`n_x = n_y`)
    and specifying **hcp** bundling will generate *rhomboidal* bundles
    (*i.e.* bundles arranged within a rhomboid) (rhombuses if
    :math:`n_x = n_y`). However, you can also enforce a specific
    *bundle geometry* which will try and reshape the bundle arrangement so
    that it "fits inside" the boundaries of a specified geometric shape.
    This allows you to generate **hcp** bundles that are trianglar,
    hexagonal, or rectangular in *shape*, as some of the examples below
    illustrate.

    To start, let's generate an hcp bundle of
    :math:`C_{\\mathrm{h}} = (10, 5)` SWCNTs and cell count
    :math:`n_x=10, n_y=3, n_z=5`.

    >>> from sknano.nanogen import NanotubeBundleGenerator
    >>> SWCNTbundle = NanotubeBundleGenerator(n=10, m=5, nx=10,
    ...                                       ny=3, nz=5)
    >>> SWCNTbundle.save_data()

    The rendered structure looks like:

    .. image:: /images/1005_hcp_10cellsx3cellsx5cells-001.png

    Now let's generate a nice hexagon bundle, 3 tubes wide, with
    :math:`C_{\\mathrm{h}} = (6, 5)`.

    >>> SWCNTbundle = NanotubeBundleGenerator(n=6, m=5, nz=5,
    ...                                       bundle_geometry='hexagon')
    >>> SWCNTbundle.save_data()

    which looks like:

    .. image:: /images/0605_hcp_7tube_hexagon-001.png

    Remember, all of the :py:meth:`~NanotubeBundleGenerator.save_data`
    methods allow you to rotate the structure data before saving:

    >>> SWCNTbundle.save_data(fname='0605_hcp_7tube_hexagon_rot-30deg.xyz',
    ...                       rot_axis='z', rotation_angle=30)

    .. image:: /images/0605_hcp_7tube_hexagon_rot-30deg-001.png

    Now, just because we can, let's make a big ass hexagon bundle with
    :math:`C_{\\mathrm{h}} = (10, 0)`.

    >>> BIGASSHEXABUN = NanotubeBundleGenerator(n=10, m=0, nx=25,
    ...                                         ny=25, nz=1,
    ...                                         bundle_geometry='hexagon')

    You're looking at 469 :math:`(10, 0)` unit cells! That's
    :math:`N_{\\mathrm{atoms}} = 18760`.

    .. image:: /images/1000_hcp_469tube_hexagon-001.png

    """

    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond, vdw_spacing=3.4,
                 bundle_packing=None, bundle_geometry=None, Lx=None, Ly=None,
                 Lz=None, fix_Lz=False, autogen=True, verbose=False):

        super(NanotubeBundleGenerator, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz, element1=element1,
            element2=element2, Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz, bond=bond,
            autogen=False, verbose=verbose)

        self.compute_bundle_params()

        self._r1 = np.zeros(3)
        self._r1[0] = \
            Nanotube.compute_dt(n=n, m=m, bond=bond, with_units=False) + \
            vdw_spacing

        self._r2 = np.zeros(3)

        if bundle_packing is None and \
                bundle_geometry in ('square', 'rectangle'):
            bundle_packing = 'cubic'
        elif bundle_packing is None:
            bundle_packing = 'hexagonal'
        elif (bundle_packing == 'cubic' and bundle_geometry not in
                (None, 'square', 'rectangle')) or \
                (bundle_packing == 'hexagonal' and bundle_geometry not in
                    (None, 'triangle', 'hexagon', 'rhombus', 'rhomboid')):
            bundle_geometry = None

        if bundle_packing == 'cubic':
            self._r2[1] = self._r1[0]
        else:
            self._r2[0] = self._r1[0] * np.cos(2 * np.pi / 3)
            self._r2[1] = self._r1[0] * np.sin(2 * np.pi / 3)

        self._bundle_packing = bundle_packing
        self._bundle_geometry = bundle_geometry

        if autogen:
            super(NanotubeBundleGenerator, self).generate_unit_cell()
            self.generate_structure_data()

    def generate_structure_data(self):
        """Generate structure data."""
        super(NanotubeBundleGenerator, self).generate_structure_data()

        self._Ntubes = 0

        swcnt0 = copy.deepcopy(self.structure_atoms)
        self.structure_atoms = Atoms()
        if self._bundle_geometry == 'hexagon':
            nrows = max(self._nx, self._ny, 3)
            if nrows % 2 != 1:
                nrows += 1

            ntubes_per_end_rows = int((nrows + 1) / 2)

            row = 0
            ntubes_per_row = nrows
            while ntubes_per_row >= ntubes_per_end_rows:
                if row == 0:
                    for n in xrange(ntubes_per_row):
                        swcnt = Atoms(atoms=swcnt0, deepcopy=True)
                        swcnt.center_CM()
                        dr = n * self._r1
                        swcnt.translate(dr)
                        self.structure_atoms.extend(swcnt.atoms)
                        self._Ntubes += 1
                else:
                    for nx in xrange(1, ntubes_per_row + 1):
                        for ny in (-row, row):
                            swcnt = Atoms(atoms=swcnt0, deepcopy=True)
                            swcnt.center_CM()
                            dy = np.zeros(3)
                            dy[0] = abs(ny) * self._r2[0]
                            dy[1] = ny * self._r2[1]
                            dr = nx * self._r1 + dy
                            swcnt.translate(dr)
                            self.structure_atoms.extend(swcnt.atoms)
                            self._Ntubes += 1
                row += 1
                ntubes_per_row = nrows - row
        else:
            for nx in xrange(self._nx):
                for ny in xrange(self._ny):
                    swcnt = Atoms(atoms=swcnt0, deepcopy=True)
                    swcnt.center_CM()
                    dr = nx * self._r1 + ny * self._r2
                    swcnt.translate(dr)
                    self.structure_atoms.extend(swcnt.atoms)
                    self._Ntubes += 1
        self._Natoms_per_bundle = \
            self.compute_Natoms_per_bundle(n=self._n, m=self._m,
                                           nz=self._nz,
                                           Ntubes=self._Ntubes)

    def save_data(self, fname=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True):
        """Save structure data.

        Parameters
        ----------
        fname : {None, str}, optional
            file name string
        structure_format : {None, str}, optional
            chemical file format of saved structure data. Must be one of:

                - xyz
                - data

            If ``None``, then guess based on ``fname`` file extension or
            default to ``xyz`` format.
        rotation_angle : {None, float}, optional
            Angle of rotation
        rot_axis : {'x', 'y', 'z'}, optional
            Rotation axis
        deg2rad : bool, optional
            Convert ``rotation_angle`` from degrees to radians.
        center_CM : bool, optional
            Center center-of-mass on origin.

        """
        if (fname is None and structure_format not in
                supported_structure_formats) or \
                (fname is not None and not
                    fname.endswith(supported_structure_formats) and
                    structure_format not in supported_structure_formats):
            structure_format = default_structure_format

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
            fname += '.' + structure_format

        else:
            if fname.endswith(supported_structure_formats) and \
                    structure_format is None:
                for ext in supported_structure_formats:
                    if fname.endswith(ext):
                        structure_format = ext
                        break
            else:
                if structure_format is None or \
                        structure_format not in supported_structure_formats:
                    structure_format = default_structure_format

        super(NanotubeBundleGenerator, self).save_data(
            fname=fname, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            deg2rad=deg2rad, center_CM=center_CM)


class MWNTGenerator(NanotubeGenerator, NanotubeBundle):
    u"""Class for generating multi-walled nanotubes.

    .. versionadded:: 0.2.7

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis atoms 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    vdw_spacing : float, optional
        van der Waals distance between nearest neighbor tubes
    bundle_packing : {None, 'hexagonal', 'cubic'}, optional
        close packing arrangement of bundles
    bundle_geometry : {None, 'triangle', 'hexagon', 'square', 'rectangle',
                       'rhombus', 'rhomboid'}, optional
    Lx, Ly, Lz : float, optional
        length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.
    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.
    max_shells : int, optional
        Maximum number of shells per MWNT.
    min_shell_diameter : float, optional
        Minimum shell diameter, in units of **Angstroms**.
    shell_spacing : float, optional
        Shell spacing in units of **Angstroms**. Default
        value is the van der Waals interaction distance of 3.4 Angstroms.
    inner_shell_Ch_type : {None, 'armchair', AC', 'zigzag', 'ZZ', 'achiral',
                           'chiral'}, optional
        If `None`, the chiralities of the inner shells are constrained only
        by their diameter and will be chosen randomly if more than one
        candidate chirality exists. If not `None`, then the inner
        shell chirality type will be added as a constraint.
    autogen : bool, optional
        if `True`, automatically call
        :py:meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :py:meth:`~NanotubeGenerator.generate_structure_data`.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.nanogen import MWNTGenerator
    >>> mwnt = MWNTGenerator(n=40, m=40, max_shells=5, Lz=1.0, fix_Lz=True)
    >>> mwnt.save_data()

    .. image:: /images/5shell_mwnt_4040_outer_Ch_1cellx1cellx4.06cells-001.png

    >>> mwntbundle = MWNTGenerator(n=40, m=40, max_shells=5, Lz=1.0,
    ...                            fix_Lz=True, bundle_geometry='hexagon')
    >>> mwntbundle.save_data()

    .. image:: /images/5shell_mwnt_4040_outer_Ch_hcp_7tube_hexagon-001.png

    """
    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond, vdw_spacing=3.4,
                 bundle_packing=None, bundle_geometry=None, Lx=None, Ly=None,
                 Lz=None, fix_Lz=False, max_shells=None,
                 min_shell_diameter=0.0, shell_spacing=3.4,
                 inner_shell_Ch_type=None, autogen=True, verbose=False):

        super(MWNTGenerator, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz, element1=element1,
            element2=element2, Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz, bond=bond,
            autogen=False, verbose=verbose)

        self._Lzmin = np.inf

        self._max_shells = max_shells
        if max_shells is None:
            self._max_shells = np.inf

        self._min_shell_diameter = min_shell_diameter
        self._shell_spacing = shell_spacing
        self._inner_shell_Ch_type = inner_shell_Ch_type

        self._Nshells_per_tube = 1
        self._Natoms_per_tube = 0

        self.compute_bundle_params()

        self._r1 = np.zeros(3)
        self._r1[0] = Nanotube.compute_dt(n=n, m=m, bond=bond) + \
            vdw_spacing

        self._r2 = np.zeros(3)

        if bundle_packing is None and \
                bundle_geometry in ('square', 'rectangle'):
            bundle_packing = 'cubic'
        elif bundle_packing is None:
            bundle_packing = 'hexagonal'
        elif (bundle_packing == 'cubic' and bundle_geometry not in
                (None, 'square', 'rectangle')) or \
                (bundle_packing == 'hexagonal' and bundle_geometry not in
                    (None, 'triangle', 'hexagon', 'rhombus', 'rhomboid')):
            bundle_geometry = None

        if bundle_packing == 'cubic':
            self._r2[1] = self._r1[0]
        else:
            self._r2[0] = self._r1[0] * np.cos(2 * np.pi / 3)
            self._r2[1] = self._r1[0] * np.sin(2 * np.pi / 3)

        self._bundle_packing = bundle_packing
        self._bundle_geometry = bundle_geometry

        if autogen:
            super(MWNTGenerator, self).generate_unit_cell()
            self.generate_structure_data()

    def generate_unit_cell(self, n=int, m=int):
        """Generate the unit cell of a MWNT shell"""
        eps = 0.01
        bond = self._bond
        e1 = self._element1
        e2 = self._element2

        N = Nanotube.compute_N(n=n, m=m)
        aCh = Nanotube.compute_chiral_angle(n=n, m=m, rad2deg=False)
        rt = Nanotube.compute_rt(n=n, m=m, bond=bond, with_units=False)
        T = Nanotube.compute_T(n=n, m=m, bond=bond, with_units=False)

        tau = Nanotube.compute_tau(n=n, m=m, bond=bond, with_units=False)
        dtau = bond * np.sin(np.pi / 6 - aCh)

        psi = Nanotube.compute_psi(n=n, m=m)
        dpsi = bond * np.cos(np.pi / 6 - aCh) / rt

        unit_cell = Atoms()

        for i in xrange(1, N + 1):
            x1 = rt * np.cos(i * psi)
            y1 = rt * np.sin(i * psi)
            z1 = i * tau

            while z1 > T - eps:
                z1 -= T

            atom1 = Atom(e1, x=x1, y=y1, z=z1)
            atom1.rezero_coords()

            unit_cell.append(atom1)

            x2 = rt * np.cos(i * psi + dpsi)
            y2 = rt * np.sin(i * psi + dpsi)
            z2 = i * tau - dtau
            while z2 > T - eps:
                z2 -= T

            atom2 = Atom(e2, x=x2, y=y2, z=z2)
            atom2.rezero_coords()

            unit_cell.append(atom2)

        return unit_cell

    def generate_structure_data(self):
        """Generate structure data."""
        super(MWNTGenerator, self).generate_structure_data()

        self._Ntubes = 0

        dt = []
        Ch = []
        for n in xrange(0, 201):
            for m in xrange(0, 201):
                dt.append(Nanotube.compute_dt(n=n, m=m, bond=self.bond))
                Ch.append((n, m))
        dt = np.asarray(dt)
        Ch = np.asarray(Ch)

        swnt0 = copy.deepcopy(self.structure_atoms)
        mwnt0 = Atoms(atoms=swnt0, deepcopy=True)
        self._Lzmin = min(self._Lzmin, self._Lz)
        mwnt0.center_CM()

        next_dt = self.dt - 2 * self._shell_spacing
        while next_dt >= self._min_shell_diameter and \
                self._Nshells_per_tube < self._max_shells:
            # get chiral indices for next_dt
            next_Ch_candidates = []
            delta_dt = 0.05
            while len(next_Ch_candidates) == 0:
                if self._inner_shell_Ch_type in ('AC', 'armchair'):
                    next_Ch_candidates = \
                        Ch[np.where(
                            np.logical_and(np.abs(dt - next_dt) <= delta_dt,
                                           Ch[:,0] == Ch[:,1]))]
                elif self._inner_shell_Ch_type in ('ZZ', 'zigzag'):
                    next_Ch_candidates = \
                        Ch[np.where(
                            np.logical_and(np.abs(dt - next_dt) <= delta_dt,
                                           np.logical_or(Ch[:,0] == 0,
                                                         Ch[:,1] == 0)))]
                elif self._inner_shell_Ch_type == 'achiral':
                    next_Ch_candidates = \
                        Ch[np.where(
                            np.logical_and(np.abs(dt - next_dt) <= delta_dt,
                                           np.logical_or(
                                               Ch[:,0] == Ch[:,1],
                                               np.logical_or(
                                                   Ch[:,0] == 0,
                                                   Ch[:,1] == 0))))]
                elif self._inner_shell_Ch_type == 'chiral':
                    next_Ch_candidates = \
                        Ch[np.where(
                            np.logical_and(np.abs(dt - next_dt) <= delta_dt,
                                           np.logical_and(
                                               Ch[:,0] != Ch[:,1],
                                               np.logical_and(
                                                   Ch[:,0] != 0,
                                                   Ch[:,1] != 1))))]
                else:
                    next_Ch_candidates = \
                        Ch[np.where(np.abs(dt - next_dt) <= delta_dt)]

                next_dt -= delta_dt

                #delta_dt += 0.05

            n, m = \
                next_Ch_candidates[
                    np.random.choice(np.arange(len(next_Ch_candidates)))]
            # generate unit cell for new chiral indices
            unit_cell = self.generate_unit_cell(n=n, m=m)
            if self._verbose:
                print('next_dt: {:.4f}'.format(next_dt))
                print('n, m = {}, {}'.format(n, m))
                print('unit_cell.Natoms: {}\n'.format(unit_cell.Natoms))
            T = Nanotube.compute_T(n=n, m=m, bond=self.bond)
            Lz = Nanotube.compute_Lz(n=n, m=m, bond=self.bond, nz=self.nz)
            self._Lzmin = min(self._Lzmin, Lz)
            shell_atoms = Atoms()
            for nz in xrange(int(np.ceil(self._nz))):
                dr = np.array([0.0, 0.0, nz * T])
                for uc_atom in unit_cell.atoms:
                    nt_atom = Atom(uc_atom.symbol)
                    nt_atom.r = uc_atom.r + dr
                    shell_atoms.append(nt_atom)
            shell_atoms.center_CM()
            mwnt0.extend(shell_atoms.atoms)
            next_dt -= 2 * self._shell_spacing
            self._Nshells_per_tube += 1

        if self._L0 is not None and self._fix_Lz:
            mwnt0.clip_bounds(abs_limit=(10 * self._L0 + 0.5) / 2,
                              r_indices=[2])
        else:
            mwnt0.clip_bounds(abs_limit=(10 * self._Lzmin + 0.5) / 2,
                              r_indices=[2])

        self._Natoms_per_tube = mwnt0.Natoms

        self.structure_atoms = Atoms()

        if self._bundle_geometry == 'hexagon':
            nrows = max(self._nx, self._ny, 3)
            if nrows % 2 != 1:
                nrows += 1

            ntubes_per_end_rows = int((nrows + 1) / 2)

            row = 0
            ntubes_per_row = nrows
            while ntubes_per_row >= ntubes_per_end_rows:
                if row == 0:
                    for n in xrange(ntubes_per_row):
                        mwnt = Atoms(atoms=mwnt0.atoms, deepcopy=True)
                        mwnt.center_CM()
                        dr = n * self._r1
                        mwnt.translate(dr)
                        self.structure_atoms.extend(mwnt.atoms)
                        self._Ntubes += 1
                else:
                    for nx in xrange(1, ntubes_per_row + 1):
                        for ny in (-row, row):
                            mwnt = Atoms(atoms=mwnt0.atoms, deepcopy=True)
                            mwnt.center_CM()
                            dy = np.zeros(3)
                            dy[0] = abs(ny) * self._r2[0]
                            dy[1] = ny * self._r2[1]
                            dr = nx * self._r1 + dy
                            mwnt.translate(dr)
                            self.structure_atoms.extend(mwnt.atoms)
                            self._Ntubes += 1
                row += 1
                ntubes_per_row = nrows - row
        else:
            for nx in xrange(self._nx):
                for ny in xrange(self._ny):
                    mwnt = Atoms(atoms=mwnt0.atoms, deepcopy=True)
                    mwnt.center_CM()
                    dr = nx * self._r1 + ny * self._r2
                    mwnt.translate(dr)
                    self.structure_atoms.extend(mwnt.atoms)
                    self._Ntubes += 1
        self._Natoms_per_bundle = self._Ntubes * self._Natoms_per_tube

        if self._verbose:
            print('Ntubes: {}'.format(self._Ntubes))
            print('Nshells_per_tube: {}'.format(self._Nshells_per_tube))
            print('Natoms_per_tube: {}'.format(self._Natoms_per_tube))
            print('Natoms_per_bundle: {}'.format(self._Natoms_per_bundle))

    def save_data(self, fname=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True):
        """Save structure data.

        Parameters
        ----------
        fname : {None, str}, optional
            file name string
        structure_format : {None, str}, optional
            chemical file format of saved structure data. Must be one of:

                - xyz
                - data

            If ``None``, then guess based on ``fname`` file extension or
            default to ``xyz`` format.
        rotation_angle : {None, float}, optional
            Angle of rotation
        rot_axis : {'x', 'y', 'z'}, optional
            Rotation axis
        deg2rad : bool, optional
            Convert ``rotation_angle`` from degrees to radians.
        center_CM : bool, optional
            Center center-of-mass on origin.

        """
        if (fname is None and structure_format not in
                supported_structure_formats) or \
                (fname is not None and not
                    fname.endswith(supported_structure_formats) and
                    structure_format not in supported_structure_formats):
            structure_format = default_structure_format

        if fname is None:

            Nshells = '{}shell_mwnt'.format(self._Nshells_per_tube)

            chirality = '{}{}_outer_Ch'.format('{}'.format(self._n).zfill(2),
                                               '{}'.format(self._m).zfill(2))
            packing = '{}cp'.format(self._bundle_packing[0])

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

                if self._nx == 1 and self._ny == 1:
                    fname_wordlist = (Nshells, chirality, cells)
                else:
                    fname_wordlist = (Nshells, chirality, packing, cells)
            else:
                fname_wordlist = \
                    (Nshells, chirality, packing, Ntube, self._bundle_geometry)

            fname = '_'.join(fname_wordlist)
            fname += '.' + structure_format

        else:
            if fname.endswith(supported_structure_formats) and \
                    structure_format is None:
                for ext in supported_structure_formats:
                    if fname.endswith(ext):
                        structure_format = ext
                        break
            else:
                if structure_format is None or \
                        structure_format not in supported_structure_formats:
                    structure_format = default_structure_format

        super(MWNTGenerator, self).save_data(
            fname=fname, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            deg2rad=deg2rad, center_CM=center_CM)
