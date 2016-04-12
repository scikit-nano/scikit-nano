# -*- coding: utf-8 -*-
"""
===============================================================================
Unrolled SWNT generator (:mod:`sknano.generators.unrolled_swnt_generator`)
===============================================================================

.. currentmodule:: sknano.generators.unrolled_swnt_generator

.. todo::

   Add methods to perform fractional translation and cartesian translation
   before structure generation.

.. todo::

   Handle different units in output coordinates.

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import copy

# import numpy as np

from sknano.core import pluralize
from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from sknano.core.crystallography import SuperCell
from sknano.core.math import Vector
from sknano.core.structures import UnrolledSWNT
from .base import NanoStructureGenerator

__all__ = ['UnrolledSWNTGenerator']


class UnrolledSWNTGenerator(NanoStructureGenerator, UnrolledSWNT):
    """Class for generating unrolled nanotube structures.

    .. versionadded:: 0.2.23

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube chiral vector
        :math:`\\mathbf{C}_{h} = n\\mathbf{a}_{1} + m\\mathbf{a}_{2} = (n, m)`.
    nx, ny, nz : int, optional
        Number of repeat unit cells in the :math:`x, y, z` dimensions
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
    Lx, Ly, Lz : float, optional
        Length of bundle in :math:`x, y, z` dimensions in **Angstroms**.
        Overrides the :math:`n_x, n_y, n_z` cell values.

        .. versionchanged:: 0.4.0

           Changed units from nanometers to **Angstroms**

    fix_Lz : bool, optional
        Generate the nanotube with length as close to the specified
        :math:`L_z` as possible. If `True`, then
        non integer :math:`n_z` cells are permitted.
    autogen : bool, optional
        if `True`, automatically call
        :meth:`~UnrolledSWNTGenerator.generate`.
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

    >>> unrolled_swnt = UnrolledSWNTGenerator(10, 5)
    >>> unrolled_swnt.save()

    The rendered structure looks like:

    .. image:: /images/unrolled_1005_1cellx1cell-2.png

    """
    def generate(self, finalize=True):
        """Generate structure data."""
        self.structure.clear()
        layer0 = Atoms()
        scaling_matrix = self.scaling_matrix
        supercell = SuperCell(self.unit_cell, scaling_matrix)

        for atom in supercell:
            layer0.append(Atom(**atom.todict()))
        layer0.center_centroid()

        lattice_shift = Vector(p0=supercell.basis.centroid,
                               p=layer0.centroid)
        lattice_shift.z = self.nlayers * lattice_shift.z
        self.lattice_shift = Vector(lattice_shift)

        for nlayer in range(self.nlayers):
            layer = copy.deepcopy(layer0)
            layer.translate(Vector([0, nlayer * self.layer_spacing, 0]))
            [setattr(atom, 'mol', nlayer + 1) for atom in layer]
            if (nlayer % 2) != 0:
                layer.translate(self.layer_shift)
            layer.rotate(angle=self.layer_rotation_angles[nlayer], axis='z')
            self.atoms.extend(layer)
            self.layers.append(layer)

        if finalize:
            self.finalize()

    @classmethod
    def generate_fname(cls, n=None, m=None, nx=None, nz=None,
                       fix_Lx=False, fix_Lz=False, **kwargs):
        fname = '{}{}'.format('{}'.format(n).zfill(2), '{}'.format(m).zfill(2))
        nx_fmtstr = '{:.2f}' if fix_Lx else '{:.0f}'
        nx = ''.join((nx_fmtstr.format(nx), pluralize('cell', nx)))
        nz_fmtstr = '{:.2f}' if fix_Lz else '{:.0f}'
        nz = ''.join((nz_fmtstr.format(nz), pluralize('cell', nz)))
        cells = 'x'.join((nx, nz))
        return '_'.join(('unrolled', fname, cells))

    def save(self, fname=None, outpath=None, structure_format=None,
             center_centroid=True, **kwargs):
        """Save structure data.

        See :meth:`~sknano.generators.GeneratorBase.save` method
        for documentation.

        """
        if fname is None:
            fname = self.generate_fname(n=self.n, m=self.m,
                                        nx=self.nx, nz=self.nz,
                                        fix_Lx=self.fix_Lx,
                                        fix_Lz=self.fix_Lz)

        if center_centroid:
            self.atoms.center_centroid()

        super().save(fname=fname, outpath=outpath,
                     structure_format=structure_format,
                     center_centroid=False, **kwargs)
