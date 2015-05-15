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

from operator import attrgetter

import numpy as np

from ._atoms import Atoms

__all__ = ['BasisAtoms']


class BasisAtoms(Atoms):
    """An `Atoms` sub-class for crystal structure basis atoms.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.BasisAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `BasisAtoms`}, optional
        if not `None`, then a list of `BasisAtom` instance objects or an
        existing `BasisAtoms` instance object.

    """
    def __init__(self, atoms=None):

        super().__init__(atoms=atoms)
