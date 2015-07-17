# -*- coding: utf-8 -*-
"""
===============================================================================
Base class for structure molecule (:mod:`sknano.core.molecules._molecule`)
===============================================================================

.. currentmodule:: sknano.core.molecules._molecule

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# from collections import OrderedDict
from functools import total_ordering

# import numbers
# import numpy as np

from sknano.core import BaseClass
# from sknano.core.math import Vector

__all__ = ['Molecule']


@total_ordering
class Molecule(BaseClass):
    """Base class for abstract representation of molecule.

    Parameters
    ----------
    atoms : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instances or an
        `Atoms` instance.

    """
    def __init__(self, atoms=None, **kwargs):
        super().__init__(**kwargs)
        self.atoms = atoms
        self.fmtstr = "atoms={atoms!r}"

    def __eq__(self, other):
        """Test equality of two `Molecule` object instances."""
        return self is other or self.atoms == other.atoms

    def __lt__(self, other):
        """Test if `self` is *less than* `other`."""
        return self.atoms < other.atoms

    def todict(self):
        return dict(atoms=self.atoms)
