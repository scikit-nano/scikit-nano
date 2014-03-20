# -*- coding: utf-8 -*-
"""
===========================================================================
Nanotube structure generators (:mod:`sknano.nanogen._nanotube_generators`)
===========================================================================

.. currentmodule:: sknano.nanogen._nanotube_generators

.. todo::

   Add methods to perform fractional translation and cartesian translation
   before structure generation.

.. todo::

   Handle different units in output coordinates.

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import copy
#import itertools
#import sys

import numpy as np

from ..chemistry import Atom, Atoms
from ..tools import plural_word_check
from ..tools.refdata import CCbond

from ._nanotubes import Nanotube
from ._structure_generator import StructureGenerator

__all__ = ['NanotubeGenerator', 'UnrolledNanotubeGenerator', 'MWNTGenerator']


class NanotubeGenerator(Nanotube, StructureGenerator):
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
        Element symbol or atomic number of basis
        :class:`~sknano.chemistry.Atoms` 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    Lz : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the `nz` value.

        .. versionadded:: 0.2.5

    tube_length : float, optional
        Length of nanotube in units of **nanometers**.
        Overrides the `nz` value.

        .. deprecated:: 0.2.5
           Use `Lz` instead

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.

        .. versionadded:: 0.2.6

    autogen : bool, optional
        if `True`, automatically call
        :meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :meth:`~NanotubeGenerator.generate_structure_data`.
    with_units : bool, optional
        Attach `units` to physical quantities
    units : None, optional
        System of units to use.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------
    First, load the :class:`~sknano.nanogen.NanotubeGenerator` class.

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
                 element1='C', element2='C', bond=CCbond, tube_length=None,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False,
                 with_units=False, units=None,
                 autogen=True, verbose=False):

        if tube_length is not None and Lz is None:
            Lz = tube_length

        super(NanotubeGenerator, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz,
            element1=element1, element2=element2,
            bond=bond, Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz,
            with_units=False, units=units, verbose=verbose)

        if autogen:
            self.generate_unit_cell()
            self.generate_structure_data()

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

        self._unit_cell = Atoms()

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

            self._unit_cell.append(atom1)

            x2 = rt * np.cos(i * psi + dpsi)
            y2 = rt * np.sin(i * psi + dpsi)
            z2 = i * tau - dtau
            while z2 > T - eps:
                z2 -= T

            atom2 = Atom(e2, x=x2, y=y2, z=z2)
            atom2.rezero_coords()

            if verbose:
                print('Basis Atom 2:\n{}'.format(atom2))

            self._unit_cell.append(atom2)

    def generate_structure_data(self):
        """Generate structure data."""
        self._structure_atoms = Atoms()
        for nz in xrange(int(np.ceil(self._nz))):
            dr = np.array([0.0, 0.0, nz * self.T])
            for uc_atom in self._unit_cell:
                nt_atom = Atom(uc_atom.symbol)
                nt_atom.r = uc_atom.r + dr
                self._structure_atoms.append(nt_atom)

    def save_data(self, fname=None, outpath=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True):
        """Save structure data.

        See :meth:`~sknano.nanogen.StructureGenerator.save_data` method
        for documentation.

        """
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

        if center_CM:
            self._structure_atoms.center_CM()

        if self._L0 is not None and self._fix_Lz:
            self._structure_atoms.clip_bounds(
                min_limit=-(10 * self._L0 + 0.2) / 2, r_indices=[2])
            if center_CM:
                self._structure_atoms.center_CM()
            self._structure_atoms.clip_bounds(
                max_limit=(10 * self._L0 + 0.2) / 2, r_indices=[2])
            if center_CM:
                self._structure_atoms.center_CM()

        super(NanotubeGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            deg2rad=deg2rad, center_CM=False)


class CappedNanotubeGenerator(Nanotube, StructureGenerator):
    u"""Class for generating capped nanotube structures."""
    pass


class UnrolledNanotubeGenerator(Nanotube, StructureGenerator):
    u"""Class for generating unrolled nanotube structures.

    .. versionadded:: 0.2.23

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.chemistry.Atoms` 1 and 2
    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    Lx, Ly, Lz : float, optional
        Length of bundle in :math:`x, y, z` dimensions in **nanometers**.
        Overrides the :math:`n_x, n_y, n_z` cell values.
    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.
    autogen : bool, optional
        if `True`, automatically call
        :meth:`~NanotubeGenerator.generate_unit_cell`,
        followed by :meth:`~NanotubeGenerator.generate_structure_data`.
    with_units : bool, optional
        Attach `units` to physical quantities
    units : None, optional
        System of units to use.
    verbose : bool, optional
        if `True`, show verbose output

    Notes
    -----
    The `UnrolledNanotubeGenerator` class generates graphene using the
    nanotube unit cell defined by the chiral vector
    :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    If you want to generate graphene with an armchair or zigzag edge using
    `length` and `width` parameters, see the
    :class:`~sknano.nanogen.GrapheneGenerator` class.

    .. seealso:: :class:`~sknano.nanogen.GrapheneGenerator`


    Examples
    --------
    First, load the :class:`~sknano.nanogen.UnrolledNanotubeGenerator`
    class.

    >>> from sknano.nanogen import UnrolledNanotubeGenerator

    Now let's generate an unrolled :math:`\\mathbf{C}_{\\mathrm{h}} = (10, 5)`
    SWCNT unit cell.

    >>> flatswcnt = UnrolledNanotubeGenerator(n=10, m=5)
    >>> flatswcnt.save_data()

    The rendered structure looks like:

    """

    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False,
                 nlayers=None, layer_spacing=None, stacking_order=None,
                 with_units=False, units=None, autogen=True, verbose=False):

        super(UnrolledNanotubeGenerator, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz,
            element1=element1, element2=element2,
            bond=bond, Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz,
            with_units=False, units=units, verbose=verbose)

        if autogen:
            self.generate_unit_cell()
            self.generate_structure_data()

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

        self._unit_cell = Atoms()

        for i in xrange(N):
            x1 = rt * i * psi
            z1 = i * tau

            while z1 > T - eps:
                z1 -= T

            atom1 = Atom(e1, x=x1, z=z1)
            atom1.rezero_coords()

            if verbose:
                print('Basis Atom 1:\n{}'.format(atom1))

            self._unit_cell.append(atom1)

            x2 = rt * (i * psi + dpsi)
            z2 = i * tau - dtau
            while z2 > T - eps:
                z2 -= T

            atom2 = Atom(e2, x=x2, z=z2)
            atom2.rezero_coords()

            if verbose:
                print('Basis Atom 2:\n{}'.format(atom2))

            self._unit_cell.append(atom2)

    def generate_structure_data(self):
        """Generate structure data."""
        self._structure_atoms = Atoms()
        for nx in xrange(self.nx):
            for nz in xrange(int(np.ceil(self.nz))):
                dr = np.array([nx * self.Ch, 0.0, nz * self.T])
                for uc_atom in self._unit_cell:
                    nt_atom = Atom(uc_atom.symbol)
                    nt_atom.r = uc_atom.r + dr
                    self._structure_atoms.append(nt_atom)

    def save_data(self, fname=None, outpath=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True):
        """Save structure data.

        See :meth:`~sknano.nanogen.StructureGenerator.save_data` method
        for documentation.

        """
        if fname is None:
            chirality = '{}{}f'.format('{}'.format(self._n).zfill(2),
                                       '{}'.format(self._m).zfill(2))
            nx = self.nx
            ny = self.ny
            fname_wordlist = None
            if nx != 1 or ny != 1:
                nx = ''.join(('{}'.format(self.nx),
                              plural_word_check('cell', self.nx)))
                ny = ''.join(('{}'.format(self.ny),
                              plural_word_check('cell', self.ny)))

                if self._assume_integer_unit_cells:
                    nz = ''.join(('{}'.format(self.nz),
                                  plural_word_check('cell', self.nz)))
                else:
                    nz = ''.join(('{:.2f}'.format(self.nz),
                                  plural_word_check('cell', self.nz)))

                cells = 'x'.join((nx, ny, nz))
                fname_wordlist = (chirality, cells)
            else:
                if self._assume_integer_unit_cells:
                    nz = ''.join(('{}'.format(self.nz),
                                  plural_word_check('cell', self.nz)))
                else:
                    nz = ''.join(('{:.2f}'.format(self.nz),
                                  plural_word_check('cell', self.nz)))

                fname_wordlist = (chirality, nz)

            fname = '_'.join(fname_wordlist)

        if center_CM:
            self._structure_atoms.center_CM()

        if self._L0 is not None and self._fix_Lz:
            self._structure_atoms.clip_bounds(
                min_limit=-(10 * self._L0 + 0.2) / 2, r_indices=[2])
            if center_CM:
                self._structure_atoms.center_CM()
            self._structure_atoms.clip_bounds(
                max_limit=(10 * self._L0 + 0.2) / 2, r_indices=[2])
            if center_CM:
                self._structure_atoms.center_CM()

        super(UnrolledNanotubeGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            deg2rad=deg2rad, center_CM=False)


class MWNTGenerator(NanotubeGenerator):
    u"""Class for generating single, multi-walled nanotubes.

    .. versionchanged:: 0.2.20

       `MWNTGenerator` no longer generates MWNT *bundles*, only *single*
       MWNTs. To generate bundled MWNT structure data, use the
       `MWNTBundleGenerator` class.

    .. versionadded:: 0.2.8

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions.
    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.chemistry.Atoms` 1 and 2
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
    autogen : bool, optional
        if `True`, automatically call
        :meth:`~MWNTGenerator.generate_unit_cell`,
        followed by :meth:`~MWNTGenerator.generate_structure_data`.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.nanogen import MWNTGenerator
    >>> mwnt = MWNTGenerator(n=40, m=40, max_shells=5, Lz=1.0, fix_Lz=True)
    >>> mwnt.save_data()

    .. image:: /images/5shell_mwnt_4040_outer_Ch_1cellx1cellx4.06cells-01.png

    """
    def __init__(self, n=int, m=int, nx=1, ny=1, nz=1,
                 element1='C', element2='C', bond=CCbond,
                 Lx=None, Ly=None, Lz=None, fix_Lz=False,
                 add_inner_shells=True, add_outer_shells=False,
                 max_shells=None, max_shell_diameter=np.inf,
                 min_shells=None, min_shell_diameter=0.0,
                 new_shell_type=None, shell_spacing=3.4,
                 with_units=False, units=None, autogen=True, verbose=False):

        super(MWNTGenerator, self).__init__(
            n=n, m=m, nx=nx, ny=ny, nz=nz,
            element1=element1, element2=element2, bond=bond,
            Lx=Lx, Ly=Ly, Lz=Lz, fix_Lz=fix_Lz,
            with_units=with_units, units=units,
            autogen=False, verbose=verbose)

        self._add_inner_shells = add_inner_shells
        self._add_outer_shells = add_outer_shells
        self._starting_shell_position = 'outer'

        self._max_shells = max_shells
        if max_shells is None:
            self._max_shells = 10
        self._max_shell_diameter = max_shell_diameter

        self._min_shells = min_shells
        if min_shells is None:
            self._min_shells = 2
        self._min_shell_diameter = min_shell_diameter

        self._new_shell_type = new_shell_type
        self._shell_spacing = shell_spacing

        self._Nshells_per_tube = 1
        self._Natoms_per_tube = 0

        if autogen:
            super(MWNTGenerator, self).generate_unit_cell()
            self.generate_structure_data()

    def _generate_unit_cell(self, n=int, m=int):
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
        """Generate structure data.

        .. todo::

           Load the diameter and chirality data from file instead of
           generating it every time.

        """
        Ch = []
        dt = []
        for n in xrange(0, 201):
            for m in xrange(0, 201):
                if (n <= 2 and m <= 2):
                    continue
                else:
                    Ch.append((n, m))
                    dt.append(Nanotube.compute_dt(n=n, m=m, bond=self._bond))
        Ch = np.asarray(Ch)
        dt = np.asarray(dt)

        super(MWNTGenerator, self).generate_structure_data()

        swnt0 = copy.deepcopy(self._structure_atoms)
        self._structure_atoms = Atoms(atoms=swnt0, deepcopy=True)
        self._structure_atoms.center_CM()

        self._max_shell_diameter = min(self._max_shell_diameter, dt.max())
        self._min_shell_diameter = max(self._min_shell_diameter, dt.min())

        Lzmin = self.Lz

        delta_dt = -2 * self._shell_spacing
        max_dt_diff = 0.05

        if self._add_outer_shells:
            delta_dt = -delta_dt
            self._starting_shell_position = 'inner'

        next_dt = self.dt + delta_dt
        while self._Nshells_per_tube < self._max_shells and \
                next_dt <= self._max_shell_diameter and \
                next_dt >= self._min_shell_diameter:

            # get chiral indices for next_dt
            next_Ch_candidates = []
            while len(next_Ch_candidates) == 0 and \
                    next_dt <= self._max_shell_diameter and \
                    next_dt >= self._min_shell_diameter:

                if self._new_shell_type in ('AC', 'armchair'):
                    next_Ch_candidates = \
                        Ch[np.where(
                            np.logical_and(np.abs(dt - next_dt) <= max_dt_diff,
                                           Ch[:,0] == Ch[:,1]))]
                elif self._new_shell_type in ('ZZ', 'zigzag'):
                    next_Ch_candidates = \
                        Ch[np.where(
                            np.logical_and(np.abs(dt - next_dt) <= max_dt_diff,
                                           np.logical_or(Ch[:,0] == 0,
                                                         Ch[:,1] == 0)))]
                elif self._new_shell_type == 'achiral':
                    next_Ch_candidates = \
                        Ch[np.where(
                            np.logical_and(np.abs(dt - next_dt) <= max_dt_diff,
                                           np.logical_or(
                                               Ch[:,0] == Ch[:,1],
                                               np.logical_or(
                                                   Ch[:,0] == 0,
                                                   Ch[:,1] == 0))))]
                elif self._new_shell_type == 'chiral':
                    next_Ch_candidates = \
                        Ch[np.where(
                            np.logical_and(np.abs(dt - next_dt) <= max_dt_diff,
                                           np.logical_and(
                                               Ch[:,0] != Ch[:,1],
                                               np.logical_and(
                                                   Ch[:,0] != 0,
                                                   Ch[:,1] != 1))))]
                else:
                    next_Ch_candidates = \
                        Ch[np.where(np.abs(dt - next_dt) <= max_dt_diff)]

                if self._add_outer_shells:
                    next_dt += max_dt_diff
                else:
                    next_dt -= max_dt_diff

            if len(next_Ch_candidates) > 0:
                n, m = next_Ch_candidates[
                    np.random.choice(np.arange(len(next_Ch_candidates)))]
                T = Nanotube.compute_T(n=n, m=m, bond=self._bond)
                Lz = Nanotube.compute_Lz(
                    n=n, m=m, bond=self._bond, nz=self._nz)
                Lzmin = min(Lzmin, Lz)

                # generate unit cell for new shell chiral indices
                shell_unit_cell = self._generate_unit_cell(n=n, m=m)

                if self._verbose:
                    print('new shell:\n'
                          'n, m = {}, {}\n'.format(n, m) +
                          'dt: {:.4f}\n'.format(next_dt) +
                          'shell_unit_cell.Natoms: {}\n'.format(
                              shell_unit_cell.Natoms))

                shell = Atoms()
                for nz in xrange(int(np.ceil(self._nz))):
                    dr = np.array([0.0, 0.0, nz * T])
                    for uc_atom in shell_unit_cell:
                        atom = Atom(uc_atom.symbol)
                        atom.r = uc_atom.r + dr
                        shell.append(atom)
                shell.center_CM()
                self._structure_atoms.extend(shell.atoms)
                self._Nshells_per_tube += 1
                next_dt += delta_dt
            else:
                break

        if self._L0 is not None and self._fix_Lz:
            self._structure_atoms.clip_bounds(
                abs_limit=(10 * self._L0 + 1) / 2, r_indices=[2])
        else:
            self._structure_atoms.clip_bounds(
                abs_limit=(10 * Lzmin + 1) / 2, r_indices=[2])

        self._Natoms_per_tube = self._structure_atoms.Natoms

        if self._verbose:
            print('Nshells_per_tube: {}'.format(self._Nshells_per_tube))
            print('Natoms_per_tube: {}'.format(self._Natoms_per_tube))

    def save_data(self, fname=None, outpath=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, deg2rad=True,
                  center_CM=True):
        """Save structure data.

        See :meth:`~sknano.nanogen.StructureGenerator.save_data` method
        for documentation.

        """
        if fname is None:
            Nshells = '{}shell_mwnt'.format(self._Nshells_per_tube)
            chirality = '{}{}_{}_Ch'.format('{}'.format(self._n).zfill(2),
                                            '{}'.format(self._m).zfill(2),
                                            self._starting_shell_position)
            fname_wordlist = None
            if self._assume_integer_unit_cells:
                nz = ''.join(('{}'.format(self._nz),
                              plural_word_check('cell', self._nz)))
            else:
                nz = ''.join(('{:.2f}'.format(self._nz),
                              plural_word_check('cell', self._nz)))
            fname_wordlist = (Nshells, chirality, nz)
            fname = '_'.join(fname_wordlist)

        super(MWNTGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            deg2rad=deg2rad, center_CM=center_CM)
