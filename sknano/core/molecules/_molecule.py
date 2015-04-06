# -*- coding: utf-8 -*-
"""
===============================================================================
Base class for structure molecule (:mod:`sknano.core.molecules._molecule`)
===============================================================================

.. currentmodule:: sknano.core.molecules._molecule

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext en'

#from collections import OrderedDict
from functools import total_ordering

#import numbers
#import numpy as np

#from sknano.core import xyz
#from sknano.core.math import Vector

__all__ = ['Molecule']


@total_ordering
class Molecule(object):
    """Base class for abstract representation of molecule.

    Parameters
    ----------
    id : int, optional
    atoms : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instances or an
        `Atoms` instance.

    """
    # private class var
    #_molecule_attrs = ['Xc', 'Yc', 'Zc']

    def __init__(self, id=0, atoms=None):
        self.id = id
        self.atoms = atoms

    def __str__(self):
        """Return nice string representation of `Molecule`."""
        return repr(self)

    def __repr__(self):
        """Return canonical string representation of `Molecule`."""
        return "Molecule(id={!r}, atoms={!r})".format(self.id, self.atoms)

    def __eq__(self, other):
        """Test equality of two `Molecule` object instances."""
        if self is other or (self.id == other.id and
                             self.atoms == other.atoms):
            return True
        return False

    def __lt__(self, other):
        """Test if `self` is *less than* `other`."""
        if self.id == other.id:
            return self.atoms < other.atoms
        return self.id < other.id
