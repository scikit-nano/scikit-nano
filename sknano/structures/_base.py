# -*- coding: utf-8 -*-
"""
==============================================================================
Base structure classes (:mod:`sknano.structures._base`)
==============================================================================

.. currentmodule:: sknano.structures._base

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

# import numbers

from sknano.core import BaseClass
from sknano.core.atoms import vdw_radius_from_basis
from sknano.core.crystallography import BaseStructure
from sknano.core.refdata import aCC, element_data

__all__ = ['NanoStructureBase']

r_CC_vdw = element_data['C']['VanDerWaalsRadius']


class NanoStructureBase(BaseStructure, BaseClass):
    """Base class for creating abstract representations of nanostructure.

    Parameters
    ----------
    basis : {:class:`~python:str`, :class:`~python:int`, list}, optional
        Element chemical symbols or atomic numbers of basis
        :class:`~sknano.core.atoms.Atom` 1 and 2

    """
    def __init__(self, *args, basis=None, bond=None, **kwargs):

        if basis is None:
            basis = ['C', 'C']
        basis = basis[:]

        if 'element1' in kwargs:
            basis[0] = kwargs['element1']
            del kwargs['element1']

        if 'element2' in kwargs:
            basis[1] = kwargs['element2']
            del kwargs['element2']

        if 'vdw_spacing' in kwargs:
            vdw_radius = kwargs['vdw_spacing'] / 2
            del kwargs['vdw_spacing']
        elif 'vdw_radius' in kwargs:
            vdw_radius = kwargs['vdw_radius']
            del kwargs['vdw_radius']
        else:
            vdw_radius = None

        if bond is None:
            bond = aCC

        super().__init__(*args, **kwargs)

        self.bond = bond
        self.basis = basis
        self.vdw_radius = vdw_radius

    @property
    def basis(self):
        """:class:`NanoStructureBase` basis atoms."""
        return self._basis

    @basis.setter
    def basis(self, value):
        self._basis = value
        try:
            [self.crystal_cell.update_basis(element, index=i, step=2) for
             i, element in enumerate(self.basis)]
        except AttributeError:
            pass

    @basis.deleter
    def basis(self):
        del self._basis

    @property
    def vdw_radius(self):
        """van der Waals radius"""
        if self._vdw_radius is not None:
            return self._vdw_radius
        else:
            return vdw_radius_from_basis(self.basis[0], self.basis[1])

    @vdw_radius.setter
    def vdw_radius(self, value):
        self._vdw_radius = value

    @property
    def vdw_distance(self):
        """van der Waals distance."""
        return 2 * self.vdw_radius

    @property
    def element1(self):
        "Basis element 1"
        return self.basis[0]

    @element1.setter
    def element1(self, value):
        self.basis[0] = value
        try:
            self.crystal_cell.update_basis(value, index=0, step=2)
        except AttributeError:
            pass

    @property
    def element2(self):
        "Basis element 2"
        return self.basis[1]

    @element2.setter
    def element2(self, value):
        self.basis[1] = value
        try:
            self.crystal_cell.update_basis(value, index=1, step=2)
        except AttributeError:
            pass
