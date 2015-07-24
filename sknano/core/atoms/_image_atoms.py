# -*- coding: utf-8 -*-
"""
===============================================================================
Container class for `ImageAtom`\ s (:mod:`sknano.core.atoms._image_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._image_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

from ._atoms import Atoms
from ._image_atom import ImageAtom

__all__ = ['ImageAtoms']


class ImageAtoms(Atoms):
    """An `Atoms` sub-class for `ImageAtom`\ s.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.ImageAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `ImageAtoms`}, optional
        if not `None`, then a list of `ImageAtom` instance objects or an
        existing `ImageAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return ImageAtom

    def sort(self, key=attrgetter('i'), reverse=False):
        super().sort(key=key, reverse=reverse)
