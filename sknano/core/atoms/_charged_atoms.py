# -*- coding: utf-8 -*-
"""
===============================================================================
Container class for `ChargedAtom`\ s (:mod:`sknano.core.atoms._charged_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._charged_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from ._atoms import Atoms
from ._charged_atom import ChargedAtom

__all__ = ['ChargedAtoms']


class ChargedAtoms(Atoms):
    """An `Atoms` sub-class for `ChargedAtom`\ s.

    A container class for `ChargedAtom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `ChargedAtoms`}, optional
        if not `None`, then a list of `ChargedAtom` instance objects or an
        existing `ChargedAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return ChargedAtom

    def sort(self, key=attrgetter('q'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def charges(self):
        """Return array of `ChargedAtom` charges."""
        return np.asarray([atom.q for atom in self])

    @property
    def q(self):
        """Return the total net charge of `ChargedAtoms`."""
        return self.charges.sum()
