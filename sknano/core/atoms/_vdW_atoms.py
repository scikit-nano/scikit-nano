# -*- coding: utf-8 -*-
"""
===============================================================================
Container class for `VanDerWaalsAtom`\ s (:mod:`sknano.core.atoms._vdW_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._vdW_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from ._atoms import Atoms
from ._vdW_atom import VanDerWaalsAtom
from sknano.core.refdata import element_data

__all__ = ['VanDerWaalsAtom', 'vdw_radius_from_basis']


class VanDerWaalsAtoms(Atoms):
    """An `Atoms` sub-class for `VanDerWaalsAtom`\ s.

    A container class for `VanDerWaalsAtom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `ChargedAtoms`}, optional
        if not `None`, then a list of `ChargedAtom` instance objects or an
        existing `ChargedAtoms` instance object.

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
