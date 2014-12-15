# -*- coding: utf-8 -*-
"""
===============================================================================
SWNT structure generator (:mod:`sknano.generators._swnt_generator`)
===============================================================================

.. currentmodule:: sknano.generators._swnt_generator

.. todo::

   Add methods to perform fractional translation and cartesian translation
   before structure generation.

.. todo::

   Handle different units in output coordinates.

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
#import six
from six.moves import range
__docformat__ = 'restructuredtext en'

import numpy as np

from sknano.core import pluralize
from sknano.core.math import Vector
from sknano.structures import SWNT
from sknano.utils.geometric_shapes import Cuboid
from ._base import Atom, Atoms, GeneratorBase

__all__ = ['SWNTGenerator']


class SWNTGenerator(SWNT, GeneratorBase):
    """Class for generating `SWNT` structures.

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
        :class:`~sknano.core.Atom` 1 and 2
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
        :meth:`~SWNTGenerator.generate_unit_cell`,
        followed by :meth:`~SWNTGenerator.generate_structure_data`.
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------
    First, load the :class:`~sknano.generators.SWNTGenerator` class.

    >>> from sknano.generators import SWNTGenerator

    Now let's generate a :math:`\\mathbf{C}_{\\mathrm{h}} = (10, 5)`
    SWCNT unit cell.

    >>> nt = SWNTGenerator(n=10, m=5)
    >>> nt.save_data(fname='10,5_unit_cell.xyz')

    The rendered structure looks like (orhographic view):

    .. image:: /images/10,5_unit_cell_orthographic_view.png

    and the perspective view:

    .. image:: /images/10,5_unit_cell_perspective_view.png

    """
    def __init__(self, autogen=True, **kwargs):

        super(SWNTGenerator, self).__init__(**kwargs)

        if autogen:
            self.generate_unit_cell()
            self.generate_structure_data()

    def generate_unit_cell(self):
        """Generate the nanotube unit cell."""
        eps = 0.01

        e1 = self.element1
        e2 = self.element2
        N = self.N
        T = self.T
        rt = self.rt

        psi, tau, dpsi, dtau = self.unit_cell_symmetry_params

        if self.verbose:
            print('dpsi: {}'.format(dpsi))
            print('dtau: {}\n'.format(dtau))

        self.unit_cell = Atoms()

        for i in range(N):
            x1 = rt * np.cos(i * psi)
            y1 = rt * np.sin(i * psi)
            z1 = i * tau

            while z1 > T - eps:
                z1 -= T

            if z1 < 0:
                z1 += T

            if self.debug:
                print('i={}: x1, y1, z1 = ({:.6f}, {:.6f}, {:.6f})'.format(
                    i, x1, y1, z1))

            atom1 = Atom(element=e1, x=x1, y=y1, z=z1)
            atom1.rezero()

            if self.verbose:
                print('Basis Atom 1:\n{}'.format(atom1))

            self.unit_cell.append(atom1)

            x2 = rt * np.cos(i * psi + dpsi)
            y2 = rt * np.sin(i * psi + dpsi)
            z2 = i * tau - dtau

            while z2 > T - eps:
                z2 -= T

            if z2 < 0:
                z2 += T

            if self.debug:
                print('i={}: x2, y2, z2 = ({:.6f}, {:.6f}, {:.6f})'.format(
                    i, x2, y2, z2))

            atom2 = Atom(element=e2, x=x2, y=y2, z=z2)
            atom2.rezero()

            if self.verbose:
                print('Basis Atom 2:\n{}'.format(atom2))

            self.unit_cell.append(atom2)

    def generate_structure_data(self):
        """Generate structure data."""
        #self.atoms = Atoms()
        self.structure_data.clear()
        for nz in range(int(np.ceil(self.nz))):
            dr = Vector([0.0, 0.0, nz * self.T])
            for uc_atom in self.unit_cell:
                nt_atom = Atom(element=uc_atom.symbol)
                nt_atom.r = uc_atom.r + dr
                self.atoms.append(nt_atom)

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
            nz = '{}' if self._assert_integer_nz else '{:.2f}'
            nz = ''.join((nz.format(self.nz), pluralize('cell', self.nz)))
            fname_wordlist = (chirality, nz)
            fname = '_'.join(fname_wordlist)

        if self.L0 is not None and self.fix_Lz:
            Lz_cutoff = 10 * self.L0 + 1
            pmin = [-np.inf, -np.inf, -Lz_cutoff]
            pmax = [np.inf, np.inf, Lz_cutoff]
            region_bounds = Cuboid(pmin=pmin, pmax=pmax)
            region_bounds.update_region_limits()

            self.atoms.clip_bounds(region_bounds, center_before_clipping=True)

        if center_CM:
            self.atoms.center_CM()

        super(SWNTGenerator, self).save_data(
            fname=fname, outpath=outpath, structure_format=structure_format,
            rotation_angle=rotation_angle, rot_axis=rot_axis,
            anchor_point=anchor_point, deg2rad=deg2rad, center_CM=False,
            savecopy=savecopy, **kwargs)
