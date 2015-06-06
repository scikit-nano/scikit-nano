# -*- coding: utf-8 -*-
"""
===============================================================================
Atoms class with a coordination number (:mod:`sknano.core.atoms._cn_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._cn_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from ._atoms import Atoms
from ._cn_atom import CNAtom

__all__ = ['CNAtoms']


class CNAtoms(Atoms):
    """An eXtended `Atoms` class.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.CNAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `CNAtoms`}, optional
        if not `None`, then a list of `CNAtom` instance objects or an
        existing `CNAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return CNAtom

    def sort(self, key=attrgetter('CN'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def coordination_numbers(self):
        """:class:`~numpy:numpy.ndarray` of `CNAtom.CN`."""
        return np.asarray([atom.CN for atom in self])
