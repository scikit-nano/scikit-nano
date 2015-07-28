# -*- coding: utf-8 -*-
"""
===============================================================================
`StructureAtoms` container class (:mod:`sknano.core.atoms._structure_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._structure_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

from ._structure_atom import StructureAtom
from ._charged_atoms import ChargedAtoms
from ._velocity_atoms import VelocityAtoms
from ._image_atoms import ImageAtoms
from ._type_atoms import TypeAtoms
from ._poav_atoms import POAVAtoms

__all__ = ['StructureAtoms']


class StructureAtoms(POAVAtoms, ChargedAtoms, VelocityAtoms, ImageAtoms,
                     TypeAtoms):
    """An `Atoms` sub-class for structure analysis.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.StructureAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `StructureAtoms`}, optional
        if not `None`, then a list of `StructureAtom` instance objects or an
        existing `StructureAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return StructureAtom

    def sort(self, key=attrgetter('element', 'Z', 'mass', 'id', 'mol', 'type',
                                  'x', 'y', 'z', 'CN'), reverse=False):
        super().sort(key=key, reverse=reverse)

    def compute_rdf(self):
        pass

    @property
    def volume(self):
        try:
            return self._volume
        except AttributeError:
            return self.bounds.volume

    @volume.setter
    def volume(self, value):
        self._volume = float(value)
