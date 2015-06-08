# -*- coding: utf-8 -*-
"""
===============================================================================
Atoms class for crystal structure basis (:mod:`sknano.core.atoms._basis_atoms`)
===============================================================================

An `Atoms` sub-class for crystal structure basis atoms.

.. currentmodule:: sknano.core.atoms._basis_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# from operator import attrgetter

import numpy as np

from ._xyz_atoms import XYZAtoms
from ._basis_atom import BasisAtom

__all__ = ['BasisAtoms']


class BasisAtoms(XYZAtoms):
    """An `Atoms` sub-class for crystal structure basis atoms.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.BasisAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `BasisAtoms`}, optional
        if not `None`, then a list of `BasisAtom` instance objects or an
        existing `BasisAtoms` instance object.

    """

    @property
    def __atom_class__(self):
        return BasisAtom

    @property
    def rs(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`Atom.r` position \
            `Vector`\ s"""
        return np.asarray([atom.rs for atom in self])

    @property
    def xs(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`x` coordinates."""
        return self.rs[:, 0]

    @property
    def ys(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`y` coordinates."""
        return self.rs[:, 1]

    @property
    def zs(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`z` coordinates."""
        return self.rs[:, 2]

    @property
    def lattice(self):
        try:
            return self[0].lattice
        except IndexError:
            return None

    @lattice.setter
    def lattice(self, value):
        [setattr(atom, 'lattice', value) for atom in self]
