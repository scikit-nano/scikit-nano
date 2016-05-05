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

from sknano.core import deprecated, deprecated_kwargs
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
    @deprecated(since='0.4.0', alternative='lattice.a', obj_type='attribute')
    def Lx(self):
        """Axis-aligned length along the `x`-axis in **Angstroms**."""
        return self.lattice.a

    @property
    @deprecated(since='0.4.0', alternative='lattice.b', obj_type='attribute')
    def Ly(self):
        """Axis-aligned length along the `y`-axis in **Angstroms**."""
        return self.lattice.b

    @property
    @deprecated(since='0.4.0', alternative='n1', obj_type='attribute')
    def nx(self):
        """Number of unit cells along the :math:`x`-axis."""
        return self.n1

    @nx.setter
    def nx(self, value):
        self.n1 = value

    @property
    @deprecated(since='0.4.0', alternative='n2', obj_type='attribute')
    def ny(self):
        """An alias for :attr:`UnrolledSWNTMixin.nlayers`."""
        return self.n2

    @ny.setter
    def ny(self, value):
        self.n2 = value

    @property
    def l1(self):
        """Axis-aligned length along the `y`-axis in **Angstroms**."""
        return self.n1 * self.Ch

    @property
    def l2(self):
        """Axis-aligned length along the `y`-axis in **Angstroms**."""
        return self.n2 * self.layer_spacing

    @property
    def n1(self):
        """Number of nanotubes along the :math:`x`-axis."""
        return self._n1

    @n1.setter
    def n1(self, value):
        """Set :math:`n_x`"""
        if not (isinstance(value, numbers.Number) or value > 0):
            raise TypeError('Expected a positive integer.')
        self._n1 = int(value)

    @property
    def n2(self):
        """Number of nanotubes along the :math:`y`-axis."""
        return self.nlayers

    @n2.setter
    def n2(self, value):
        self.nlayers = value

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

    @deprecated_kwargs(dict(nx='n1', Lx=None, fix_Lx=None), since='0.4.0')
    def __init__(self, *args, n1=1, **kwargs):
        super().__init__(*args, **kwargs)

        self.n1 = n1

        if self.nlayers > 1 and self.stacking_order == 'AB':
            chiral_angle = compute_chiral_angle(self.n, self.m, degrees=False)
            self.layer_shift.x = self.bond * np.cos(np.pi/6 - chiral_angle)
            self.layer_shift.z = -self.bond * np.sin(np.pi/6 - chiral_angle)

        self.__generate_unit_cell()
        self.scaling_matrix = \
            [int(np.ceil(self.n1)), self.nlayers, int(np.ceil(self.n3))]
        self.fmtstr = ", ".join((super().fmtstr, "n1={n1!r}"))

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

    __generate_unit_cell = generate_unit_cell

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attr_dict = super().todict()
        attr_dict.update(dict(n1=self.n1))
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
    n1 : :class:`python:int`, optional
        Number of repeat unit cells in the :math:`x` direction, along
        the *unrolled* chiral vector.
    n3 : :class:`python:int`, optional
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
    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.core.structures import UnrolledSWNT
    >>> unrolled_swnt = UnrolledSWNT(10, 5)
    >>> unrolled_swnt
    UnrolledSWNT((10, 5), n1=1, n3=1, bond=1.42, basis=['C', 'C'], nlayers=1,
    layer_spacing=3.4, stacking_order='AB')

    """
    pass
