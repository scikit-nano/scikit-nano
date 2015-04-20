# -*- coding: utf-8 -*-
"""
===============================================================================
Base crystallography classes (:mod:`sknano.core.crystallography._base`)
===============================================================================

.. currentmodule:: sknano.core.crystallography._base

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
from builtins import object
from future.utils import with_metaclass

__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractproperty
#from sknano.core.math import Points, Vectors

#import numpy as np

__all__ = ['CrystalLattice', 'CrystalStructure']


class CrystalLattice(with_metaclass(ABCMeta, object)):
    """Abstract base class for crystal lattice systems."""

    def __init__(self):
        self.pstr = ""

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.pstr.format(**self.pdict()))

    @abstractproperty
    def primitive_vectors(self):
        """Primitive lattice vectors."""
        raise NotImplementedError


class CrystalStructure(with_metaclass(ABCMeta, object)):
    """Abstract base class for crystal structures."""

    def __init__(self):
        self.pstr = ""

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.pstr.format(**self.pdict()))

    @abstractproperty
    def basis(self):
        """Crystal structure basis."""
        raise NotImplementedError
