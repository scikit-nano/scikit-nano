# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes with van der Waals radius (:mod:`sknano.core.atoms.vdW_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms.vdW_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from operator import attrgetter
import numbers
import numpy as np

from sknano.core.refdata import element_data
from .atoms import Atom, Atoms

__all__ = ['VanDerWaalsAtom', 'VanDerWaalsAtoms', 'vdw_radius_from_basis']


class VanDerWaalsAtom(Atom):
    """An `Atom` class with a van der Waals radius attribute.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    r_vdw : float, optional

    """
    def __init__(self, *args, r_vdw=0., **kwargs):

        super().__init__(*args, **kwargs)

        self.r_vdw = r_vdw
        self.fmtstr = super().fmtstr + ", r_vdw={r_vdw!r}"

    @property
    def __atoms_class__(self):
        return VanDerWaalsAtoms

    # def __eq__(self, other):
    #     return np.allclose(self.r_vdw, other.r_vdw) and super().__eq__(other)

    # def __lt__(self, other):
    #     return (self.r_vdw < other.r_vdw and super().__le__(other)) or \
    #         (self.r_vdw <= other.r_vdw and super().__lt__(other))

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return np.allclose(self.r_vdw, other.r_vdw) and super().__eq__(other)

    def __le__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.r_vdw > other.r_vdw or not super().__le__(other):
            return False
        return True

    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.r_vdw >= other.r_vdw or not super().__lt__(other):
            return False
        return True

    def __ge__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.r_vdw < other.r_vdw or not super().__ge__(other):
            return False
        return True

    def __gt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.r_vdw <= other.r_vdw or not super().__gt__(other):
            return False
        return True

    def __dir__(self):
        attrs = super().__dir__()
        attrs.append('r_vdw')
        return attrs

    @property
    def r_vdw(self):
        """Van der Waals radius :math:`r_{\\mathrm{vdW}}` in units of \
            :math:`\\AA`.
        """
        return self._r_vdw

    @r_vdw.setter
    def r_vdw(self, value):
        """Set `VanDerWaalsAtom` radius :math:`r_{\\mathrm{vdw}}`.

        Parameters
        ----------
        value : float
            van der Waals radius in units of :math:`\\AA`.

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._r_vdw = float(value)

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        super_dict = super().todict()
        super_dict.update(dict(r_vdw=self.r_vdw))
        return super_dict


class VanDerWaalsAtoms(Atoms):
    """An `Atoms` sub-class for `VanDerWaalsAtom`\ s.

    A container class for `VanDerWaalsAtom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `VanDerWaalsAtoms`}, optional
        if not `None`, then a list of `VanDerWaalsAtom` instance objects or an
        existing `VanDerWaalsAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return VanDerWaalsAtom

    def sort(self, key=attrgetter('r_vdw'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def r_vdw(self):
        """Return array of `VanDerWaalsAtom` van der Waals radii."""
        return np.asarray([atom.r_vdw for atom in self])


def vdw_radius_from_basis(*args):
    """Helper function to return the maximum van der Waals radius from a \
        tuple of atoms."""
    r_vdw = 0
    for atom in args:
        try:
            element = atom.element
        except AttributeError:
            element = atom

        try:
            r_vdw = max(r_vdw, element_data[element]['VanDerWaalsRadius'])
        except TypeError:
            try:
                r_vdw = max(r_vdw, element_data[element]['AtomicRadius'])
            except TypeError:
                pass

    return r_vdw
