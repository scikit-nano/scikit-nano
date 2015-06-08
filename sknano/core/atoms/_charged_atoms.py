# -*- coding: utf-8 -*-
"""
===============================================================================
Extended Atoms class feature set (:mod:`sknano.core.atoms._charged_atoms`)
===============================================================================

An "eXtended" `Atoms` class for structure analysis.

.. currentmodule:: sknano.core.atoms._charged_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
from builtins import zip
__docformat__ = 'restructuredtext en'

from collections import OrderedDict
from operator import attrgetter

import numpy as np

from sknano.core import dedupe, xyz
from sknano.core.math import Vector, transformation_matrix
from sknano.utils.geometric_shapes import Cuboid  # , Rectangle
from ._atoms import Atoms
from ._charged_atom import ChargedAtom

__all__ = ['ChargedAtoms']


class ChargedAtoms(Atoms):
    """An eXtended `Atoms` class.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.ChargedAtom` instances.

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
