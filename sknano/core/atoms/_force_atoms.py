# -*- coding: utf-8 -*-
"""
===============================================================================
Container class for `ForceAtom`\ s (:mod:`sknano.core.atoms._force_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._force_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

from ._atoms import Atoms
from ._force_atom import ForceAtom

__all__ = ['ForceAtoms']


class ForceAtoms(Atoms):
    """An `Atoms` sub-class for `ForceAtom`\ s.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.ForceAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `ForceAtoms`}, optional
        if not `None`, then a list of `ForceAtom` instance objects or an
        existing `ForceAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return ForceAtom

    def sort(self, key=attrgetter('f'), reverse=False):
        super().sort(key=key, reverse=reverse)
