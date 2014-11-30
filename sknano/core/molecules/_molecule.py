# -*- coding: utf-8 -*-
"""
===============================================================================
Base class for structure molecule (:mod:`sknano.core.molecules._molecule`)
===============================================================================

.. currentmodule:: sknano.core.molecules._molecule

"""
from __future__ import absolute_import, division, print_function
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
    moleculeID : int, optional
    atoms : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instances or an
        `Atoms` instance.

    """
    # private class var
    #_molecule_attrs = ['Xc', 'Yc', 'Zc']

    def __init__(self, moleculeID=0, atoms=None):
        self.moleculeID = moleculeID
        self.atoms = atoms

    def __str__(self):
        """Return nice string representation of `Molecule`."""
        return repr(self)

    def __repr__(self):
        """Return canonical string representation of `Molecule`."""
        return "Molecule(moleculeID={!r}, atoms={!r})".format(
            self.moleculeID, self.atoms)

    def __eq__(self, other):
        """Test equality of two `Molecule` object instances."""
        if self is other or (self.moleculeID == other.moleculeID and
                             self.atoms == other.atoms):
            return True
        return False

    def __lt__(self, other):
        """Test if `self` is *less than* `other`."""
        if self.moleculeID == other.moleculeID:
            return self.atoms < other.atoms
        return self.moleculeID < other.moleculeID
