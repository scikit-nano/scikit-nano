# -*- coding: utf-8 -*-
"""
==============================================================================
SWNT bundle structure class (:mod:`sknano.structures._swnt_bundle`)
==============================================================================

.. currentmodule:: sknano.structures._swnt_bundle

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._compute_funcs import compute_bundle_density
from ._nanotube_bundle import NanotubeBundleBase
from ._swnt import SWNT

__all__ = ['SWNTBundle']


class SWNTBundle(NanotubeBundleBase, SWNT):
    """SWNT bundle structure class.

    Parameters
    ----------
    *Ch : {:class:`python:tuple` or :class:`python:int`\ s}
        Either a 2-tuple of integers (i.e., *Ch = ((n, m)) or
        2 integers (i.e., *Ch = (n, m) specifying the chiral indices
        of the nanotube chiral vector
        :math:`\\mathbf{C}_h = n\\mathbf{a}_1 + m\\mathbf{a}_2 = (n, m)`.
    nx : :class:`python:int`, optional
        Number of nanotubes along the :math:`x` axis
    ny : :class:`python:int`, optional
        Number of nanotubes along the :math:`y` axis
    nz : :class:`python:int`, optional
        Number of repeat unit cells in the :math:`z` direction, along
        the *length* of the nanotube.
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

    verbose : bool, optional
        if `True`, show verbose output

    Examples
    --------

    >>> from sknano.structures import SWNT

    Create a SWNT with :math:`\\mathbf{C}_{h} = (10, 10)` chirality.

    """

    @property
    def bundle_density(self):
        return compute_bundle_density(self.n, self.m, d_vdw=self.vdw_spacing,
                                      bond=self.bond, basis=self.basis)
