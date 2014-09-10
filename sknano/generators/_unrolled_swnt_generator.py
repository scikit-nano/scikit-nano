# -*- coding: utf-8 -*-
"""
===============================================================================
Unrolled SWNT generator (:mod:`sknano.generators._unrolled_swnt_generator`)
===============================================================================

.. currentmodule:: sknano.generators._unrolled_swnt_generator

.. todo::

   Add methods to perform fractional translation and cartesian translation
   before structure generation.

.. todo::

   Handle different units in output coordinates.

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core import pluralize
from sknano.core.math import Vector
from sknano.structures import UnrolledSWNT, compute_chiral_angle
#from sknano.utils.geometric_shapes import Cuboid
from ._base import Atom, Atoms, GeneratorMixin

__all__ = ['UnrolledSWNTGenerator']


class UnrolledSWNTGenerator(UnrolledSWNT, GeneratorMixin):
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
        :class:`~sknano.core.Atom` 1 and 2
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
    verbose : bool, optional
        if `True`, show verbose output

    Notes
    -----
    The `UnrolledSWNTGenerator` class generates graphene using the
    nanotube unit cell defined by the chiral vector
    :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    If you want to generate graphene with an armchair or zigzag edge using
    `length` and `width` parameters, see the
    :class:`~sknano.generators.GrapheneGenerator` class.

    .. seealso:: :class:`~sknano.generators.GrapheneGenerator`


    Examples
    --------
    First, load the :class:`~sknano.generators.UnrolledSWNTGenerator`
    class.

    >>> from sknano.generators import UnrolledSWNTGenerator

    Now let's generate an unrolled :math:`\\mathbf{C}_{\\mathrm{h}} = (10, 5)`
    SWCNT unit cell.

    >>> flatswcnt = UnrolledSWNTGenerator(n=10, m=5)
    >>> flatswcnt.save_data()

    The rendered structure looks like:

    """

    def __init__(self, autogen=True, **kwargs):

        super(UnrolledSWNTGenerator, self).__init__(**kwargs)

        if autogen:
            self.generate_unit_cell()
            self.generate_structure_data()

    def generate_unit_cell(self):
        """Generate the nanotube unit cell."""
        eps = 0.01
        n = self.n
        m = self.m
        bond = self.bond
        M = self.M
        T = self.T
        N = self.N
        rt = self.rt
        e1 = self.element1
        e2 = self.element2
        verbose = self.verbose

        aCh = compute_chiral_angle(n=n, m=m, rad2deg=False)

        tau = M * T / N
        dtau = bond * np.sin(np.pi / 6 - aCh)

        psi = 2 * np.pi / N
        dpsi = bond * np.cos(np.pi / 6 - aCh) / rt

        if verbose:
            print('dpsi: {}'.format(dpsi))
            print('dtau: {}\n'.format(dtau))

        self.unit_cell = Atoms()

        for i in xrange(N):
            x1 = rt * i * psi
            z1 = i * tau

            while z1 > T - eps:
                z1 -= T

            atom1 = Atom(element=e1, x=x1, z=z1)
            atom1.rezero()

            if verbose:
                print('Basis Atom 1:\n{}'.format(atom1))

            self.unit_cell.append(atom1)

            x2 = rt * (i * psi + dpsi)
            z2 = i * tau - dtau
            while z2 > T - eps:
                z2 -= T

            atom2 = Atom(element=e2, x=x2, z=z2)
            atom2.rezero()

            if verbose:
                print('Basis Atom 2:\n{}'.format(atom2))

            self.unit_cell.append(atom2)

    def generate_structure_data(self):
        """Generate structure data."""
        self.atoms = Atoms()
        for nx in xrange(self.nx):
            for nz in xrange(int(np.ceil(self.nz))):
                dr = Vector([nx * self.Ch, 0.0, nz * self.T])
                for uc_atom in self.unit_cell:
                    nt_atom = Atom(element=uc_atom.symbol)
                    nt_atom.r = uc_atom.r + dr
                    self.atoms.append(nt_atom)

    def save_data(self, fname=None, outpath=None, structure_format=None,
                  rotation_angle=None, rot_axis=None, anchor_point=None,
                  deg2rad=True, center_CM=True, savecopy=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorMixin.save_data` method
        for documentation.

        """
        if fname is None:
            chirality = '{}{}'.format('{}'.format(self.n).zfill(2),
                                      '{}'.format(self.m).zfill(2))
            nx = self.nx
            ny = self.ny
            fname_wordlist = None
            if nx != 1 or ny != 1:
                nx = ''.join(('{}'.format(self.nx),
                              pluralize('cell', self.nx)))
                ny = ''.join(('{}'.format(self.ny),
                              pluralize('cell', self.ny)))

                if self._assume_integer_unit_cells:
                    nz = ''.join(('{}'.format(self.nz),
                                  pluralize('cell', self.nz)))
                else:
                    nz = ''.join(('{:.2f}'.format(self.nz),
                                  pluralize('cell', self.nz)))

                cells = 'x'.join((nx, ny, nz))
                fname_wordlist = (chirality, cells)
            else:
                if self._assume_integer_unit_cells:
                    nz = ''.join(('{}'.format(self.nz),
                                  pluralize('cell', self.nz)))
                else:
                    nz = ''.join(('{:.2f}'.format(self.nz),
                                  pluralize('cell', self.nz)))

                fname_wordlist = (chirality, nz)

            fname = 'unrolled_' + '_'.join(fname_wordlist)

        super(UnrolledSWNTGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            anchor_point=anchor_point, deg2rad=deg2rad, center_CM=False,
            savecopy=savecopy, **kwargs)
