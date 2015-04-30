# -*- coding: utf-8 -*-
"""
===============================================================================
Base crystallography classes (:mod:`sknano.core.crystallography._base`)
===============================================================================

.. currentmodule:: sknano.core.crystallography._base

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
# from builtins import object
from builtins import super
# from future.utils import with_metaclass

__docformat__ = 'restructuredtext en'

# from abc import ABCMeta, abstractproperty
from enum import Enum

from sknano.core.math import Points, Vector, Vectors

# import numpy as np

__all__ = ['CrystalLattice', 'CrystalStructure']


class ReciprocalLatticeVectors:
    pass


class ReciprocalLattice:
    pass


class CrystalCell:
    pass


class CrystalLattice:
    """Base class for crystal lattice systems."""

    def __init__(self, a=None, b=None, c=None,
                 alpha=None, beta=None, gamma=None,
                 a1=None, a2=None, a3=None, cell_matrix=None):

        self.cell_matrix = None
        self.ortho_matrix = None
        self.orientation_matrix = None
        self.fractional_matrix = None
        self.lattice_type = None

        self.a1 = Vector()
        self.a2 = Vector()
        self.a3 = Vector()
        # self.generate_cell_vectors()


        if None not in (a1, a2, a3):
            pass
        elif cell_matrix is not None:
            pass
        else:
            self.a = a
            self.b = b
            self.c = c
            self.alpha = alpha
            self.beta = beta
            self.gamma = gamma

        self.pstr = "a={a!r}, b={b!r}, c={c!r}, " + \
            "alpha={alpha!r}, beta={beta!r}, gamma={gamma!r}"

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.pstr.format(**self.pdict()))

    def pdict(self):
        """Return `dict` of `CrystalLattice` parameters."""
        return dict(a=self.a, b=self.b, c=self.c,
                    alpha=self.alpha, beta=self.beta, gamma=self.gamma)

    # @property
    # def α(self):
    #     return self.alpha

    @property
    def alpha(self):
        """Angle between lattice vectors :math:`\\mathbf{b}` and \
        :math:`\\mathbf{c}`."""
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = value

    @alpha.deleter
    def alpha(self):
        del self._alpha

    @property
    def beta(self):
        """Angle between lattice vectors :math:`\\mathbf{c}` and \
        :math:`\\mathbf{a}`."""
        return self._beta

    @beta.setter
    def beta(self, value):
        self._beta = value

    @beta.deleter
    def beta(self):
        del self._beta

    @property
    def gamma(self):
        """Angle between lattice vectors :math:`\\mathbf{a}` and \
        :math:`\\mathbf{b}`."""
        return self._gamma

    @gamma.setter
    def gamma(self, value):
        self._gamma = value

    @gamma.deleter
    def gamma(self):
        del self._deleter

    @property
    def unit_cell(self):
        pass

    def fractional_to_cartesian(self):
        pass

    def cartesian_to_fractional(self):
        pass

    def wrap_fractional_coordinate(self):
        pass

    def wrap_cartesian_coordinate(self):
        pass

    def offset(self):
        pass

    def generate_cell_vectors(self):
        pass
        # self.a1.x = self.a
        # self.a2.x = self.b * np.cos(self.gamma)

    def space_group(self):
        pass

    def lattice_type(self):
        pass

    def lattice_vectors(self):
        pass

    def cell_vectors(self):
        pass

    def cell_matrix(self):
        pass

    def ortho_matrix(self):
        pass

    def orientation_matrix(self):
        pass

    def fractional_matrix(self):
        pass

    def cell_volume(self):
        pass

    def origin(self):
        pass


class CrystalStructure(CrystalLattice):
    """Abstract base class for crystal structures."""

    def __init__(self, basis, **kwargs):
        super().__init__(**kwargs)
        self.basis = basis

        self.pstr = "basis={basis!r}"
        for k, v in kwargs.items():
            self.pstr += ", {}={{{key!s}!r}}".format(k, key=k)

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.pstr.format(**self.pdict()))

    @property
    def basis(self):
        """Crystal structure basis."""
        return self._basis

    @basis.setter
    def basis(self, value):
        self._basis = value

    @basis.deleter
    def basis(self):
        del self._basis

    @property
    def unit_cell(self):
        pass

    def pdict(self):
        """Return `dict` of `CrystalStructure` parameters."""
        super_pdict = super().pdict()
        super_pdict.update(dict(basis=self.basis))
        return super_pdict
