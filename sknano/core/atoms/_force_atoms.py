# -*- coding: utf-8 -*-
"""
===============================================================================
Extended Atoms class feature set (:mod:`sknano.core.atoms._force_atoms`)
===============================================================================

An "eXtended" `Atoms` class for structure analysis.

.. currentmodule:: sknano.core.atoms._force_atoms

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
from ._force_atom import ForceAtom

__all__ = ['ForceAtoms']


class ForceAtoms(Atoms):
    """An eXtended `Atoms` class.

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
        super().sort(reverse=reverse)
