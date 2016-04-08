# -*- coding: utf-8 -*-
"""
==============================================================================
Unrolled SWNT structure class (:mod:`sknano.core.structures.unrolled_swnt`)
==============================================================================

.. currentmodule:: sknano.core.structures.unrolled_swnt

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numbers

import numpy as np

from sknano.core.atoms import BasisAtom as Atom, BasisAtoms as Atoms
from sknano.core.crystallography import Crystal3DLattice, UnitCell
# from sknano.core.math import Vector
# from sknano.core.refdata import aCC

from .base import NanoStructureBase
from .swnt import SWNTBase, compute_Ch, compute_T, compute_chiral_angle
from .graphene import GrapheneBase


__all__ = ['UnrolledSWNTMixin', 'UnrolledSWNTBase', 'UnrolledSWNT']


class UnrolledSWNTMixin:
    """Mixin class for unrolled nanotubes."""

    @property
    def fix_Lx(self):
        """:class:`~python:bool` indicating whether \
            :attr:`UnrolledSWNTMixin.Lx` is fixed or calculated."""
        return self._fix_Lx

    @fix_Lx.setter
    def fix_Lx(self, value):
        if not isinstance(value, bool):
            raise TypeError('Expected `True` or `False`')
        self._fix_Lx = value
        self._integral_nx = False if self.fix_Lx else True

    @property
    def Lx(self):
        """Axis-aligned length along the `x`-axis in **Angstroms**."""
        return self.nx * self.Ch

    @property
    def nx(self):
        """Number of unit cells along the :math:`x`-axis."""
        return self._nx

    @nx.setter
    def nx(self, value):
        """Set :math:`n_x`"""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a real positive number.')
        self._nx = int(value)

    @property
    def Ly(self):
        """Axis-aligned length along the `y`-axis in **Angstroms**."""
        return self.ny * self.layer_spacing

    @property
    def ny(self):
        """An alias for :attr:`UnrolledSWNTMixin.nlayers`."""
        return self._nlayers

    @ny.setter
    def ny(self, value):
        self.nalyers = value

    @property
    def nlayers(self):
        """Number of layers."""
        return self._nlayers

    @nlayers.setter
    def nlayers(self, value):
        """Set :attr:`nlayers`."""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a real positive number.')
        self._nlayers = int(value)


class UnrolledSWNTBase(UnrolledSWNTMixin, SWNTBase, GrapheneBase,
                       NanoStructureBase):
    """Base unrolled SWNT structure class."""

    def __init__(self, *args, nx=1, Lx=None, fix_Lx=False, **kwargs):

        super().__init__(*args, **kwargs)

        self.fix_Lx = fix_Lx
        if Lx is not None:
            self.nx = float(Lx) / self.Ch
        elif nx is not None:
            self.nx = nx
        else:
            self.nx = 1

        if self.nlayers > 1 and self.stacking_order == 'AB':
            chiral_angle = compute_chiral_angle(self.n, self.m, degrees=False)
            self.layer_shift.x = self.bond * np.cos(np.pi/6 - chiral_angle)
            self.layer_shift.z = -self.bond * np.sin(np.pi/6 - chiral_angle)

        self.generate_unit_cell()
        self.scaling_matrix = \
            [int(np.ceil(self.nx)), self.ny, int(np.ceil(self.nz))]
        self.fmtstr = ", ".join((super().fmtstr, "nx={nx!r}", "Lx={Lx!r}"))

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attr_dict = super().todict()
        attr_dict.update(dict(nx=self.nx, Lx=self.Lx, fix_Lx=self.fix_Lx))
        return attr_dict


class UnrolledSWNT(UnrolledSWNTBase, NanoStructureBase):
    """Unrolled SWNT structure class.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of integers (i.e., *Ch = ((n, m)) or
        2 integers (i.e., *Ch = (n, m) specifying the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nx : :class:`python:int`, optional
        Number of repeat unit cells in the :math:`x` direction, along
        the *unrolled* chiral vector.
    nz : :class:`python:int`, optional
        Number of repeat unit cells in the :math:`z` direction, along
        the *length* of the nanotube.
    basis : {:class:`python:list`}, optional
        List of :class:`python:str`\ s of element symbols or atomic number
        of the two atom basis (default: ['C', 'C'])

        .. versionadded:: 0.3.10

    element1, element2 : {str, int}, optional
        Element symbol or atomic number of basis
        :class:`~sknano.core.Atom`\ s 1 and 2

        .. deprecated:: 0.3.10
           Use `basis` instead

    bond : float, optional
        :math:`\\mathrm{a}_{\\mathrm{CC}} =` distance between
        nearest neighbor atoms. Must be in units of **Angstroms**.
    nlayers : int, optional
        Number of layers (default: 1)
    layer_spacing : float, optional
        Distance between layers in **Angstroms** (default: 3.4).
    stacking_order : {'AA', 'AB'}, optional
        Stacking order of layers.
    layer_rotation_angles : list, optional
        list of rotation angles for each layer in **degrees** if
        `degrees` is `True` (default), otherwise in radians.
        The list length must equal the number of layers.
    layer_rotation_increment : float, optional
        incremental layer rotation angle in **degrees** if
        `degrees` is `True` (default), otherwise in radians.
        Each subsequent layer will
        be rotated by `layer_rotation_increment` relative to the layer
        below it.
    Lx : float, optional
        Length of the unrolled swnt sheet along the chiral vector
        in units of **Angstroms**. Overrides the `nx` value.

        .. versionchanged:: 0.4.0

           Changed units from nanometers to **Angstroms**

    fix_Lx : bool, optional
        Generate the unrolled swnt sheet with the length along the
        chiral vector as close to the specified :math:`L_x` as possible.
        If `True`, then non integer :math:`n_x` cells are permitted.
    Lz : float, optional
        Length of the unrolled swnt sheet along the translation vector
        in units of **Angstroms**. Overrides the `nz` value.

        .. versionchanged:: 0.4.0

           Changed units from nanometers to **Angstroms**

    fix_Lz : bool, optional
        Generate the unrolled swnt sheet with the length along the
        translation vector as close to the specified :math:`L_z` as possible.
        If `True`, then non integer :math:`n_z` cells are permitted.

    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.core.structures import UnrolledSWNT
    >>> unrolled_swnt = UnrolledSWNT(10, 5)
    >>> unrolled_swnt
    UnrolledSWNT((10, 5), nx=1, nz=1, bond=1.42, basis=['C', 'C'], nlayers=1,
    layer_spacing=3.4, stacking_order='AB')

    """
    def generate_unit_cell(self):
        """Generate the nanotube unit cell."""
        eps = 0.01

        e1 = self.element1
        e2 = self.element2
        N = self.N
        T = self.T
        rt = self.rt

        psi, tau, dpsi, dtau = self.unit_cell_symmetry_params

        a = compute_Ch(self.n, self.m, bond=self.bond)
        b = self.layer_spacing
        c = compute_T(self.n, self.m, bond=self.bond, length=True)
        lattice = Crystal3DLattice.orthorhombic(a, b, c)

        basis = Atoms()
        if self.verbose:
            print('dpsi: {}'.format(dpsi))
            print('dtau: {}\n'.format(dtau))

        for i in range(N):
            for j, element in enumerate((e1, e2), start=1):
                theta = i * psi
                h = i * tau

                if j == 2:
                    theta += dpsi
                    h -= dtau

                x = rt * theta
                z = h

                while z > T - eps:
                    z -= T

                if z < 0:
                    z += T

                xs, ys, zs = \
                    lattice.cartesian_to_fractional([x, 0, z])
                if self.wrap_coords:
                    xs, ys, zs = \
                        lattice.wrap_fractional_coordinate([xs, ys, zs])

                if self.debug:
                    print('i={}: x, z = ({:.6f}, {:.6f})'.format(i, x, z))

                atom = Atom(element, lattice=lattice, xs=xs, ys=ys, zs=zs)
                atom.rezero()

                if self.verbose:
                    print('Basis Atom:\n{}'.format(atom))

                basis.append(atom)

        self.unit_cell = UnitCell(lattice=lattice, basis=basis)
