# -*- coding: utf-8 -*-
"""
==============================================================================
POAV analysis (:mod:`sknano.utils.analysis._poav_analysis`)
==============================================================================

.. currentmodule:: sknano.utils.analysis._poav_analysis

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

__all__ = ['POAV', 'POAV1', 'POAV2', 'POAVR']

import functools
import operator

import numpy as np

from sknano.core.math import vector as vec


class POAV(object):
    """Base class for POAV analysis.

    Parameters
    ----------
    sigma_bonds : `~sknano.core.atoms.Bonds`
        `~sknano.core.atoms.Bonds` instance.

    """
    def __init__(self, sigma_bonds):
        self.bonds = sigma_bonds
        self.b1 = self.bonds[0].vector
        self.b2 = self.bonds[1].vector
        self.b3 = self.bonds[2].vector

        self._v1 = self.b1
        self._v2 = self.b2
        self._v3 = self.b3

        self._T = vec.scalar_triple_product(self.V1, self.V2, self.V3) / 6

        self._pyramidalization_angles = None
        self._sigma_pi_angles = None
        self._misalignment_angles = None

    @property
    def v1(self):
        return self._v1

    @property
    def v2(self):
        return self._v2

    @property
    def v3(self):
        return self._v3

    @property
    def t(self):
        return vec.scalar_triple_product(self.v1, self.v2, self.v3)

    @property
    def vpi(self):
        return (self.reciprocal_v1 +
                self.reciprocal_v2 +
                self.reciprocal_v3) / self.t

    @property
    def Vpi(self):
        return self.vpi.unit_vector

    @property
    def reciprocal_v1(self):
        return vec.cross(self.v2, self.v3)

    @property
    def reciprocal_v2(self):
        return vec.cross(self.v3, self.v1)

    @property
    def reciprocal_v3(self):
        return vec.cross(self.v1, self.v2)

    @property
    def V1(self):
        return self.b1.unit_vector

    @property
    def V2(self):
        return self.b2.unit_vector

    @property
    def V3(self):
        return self.b3.unit_vector

    @property
    def R1(self):
        return self.b1.length

    @property
    def R2(self):
        return self.b2.length

    @property
    def R3(self):
        return self.b3.length

    @property
    def T(self):
        return self._T

    @property
    def A(self):
        return self.vpi.magnitude

    @property
    def H(self):
        return 3 * self.T / self.A

    @property
    def sigma_pi_angles(self):
        return self._sigma_pi_angles

    @sigma_pi_angles.setter
    def sigma_pi_angles(self, value):
        if not isinstance(value, list):
            raise TypeError('Expected a number')
        self._sigma_pi_angles = value

    @property
    def pyramidalization_angles(self):
        """Return the pyramidalization angle :math:`\\theta_P`"""
        return self._pyramidalization_angles

    @pyramidalization_angles.setter
    def pyramidalization_angles(self, value):
        if not isinstance(value, list):
            raise TypeError('Expected a number')
        self._pyramidalization_angles = value

    @property
    def misalignment_angles(self):
        return self._misalignment_angles

    @misalignment_angles.setter
    def misalignment_angles(self, value):
        if not isinstance(value, list):
            raise TypeError('Expected a number')
        self._misalignment_angles = value


class POAV1(POAV):
    """:class:`POAV` sub-class for POAV1 analysis."""

    def __init__(self, *args):
        super(POAV1, self).__init__(*args)

        self._v1 = self.V1
        self._v2 = self.V2
        self._v3 = self.V3

    @property
    def m(self):
        cos2sigmapi = np.cos(np.mean(self.sigma_pi_angles))**2
        return 2 * cos2sigmapi / (1 - 3 * cos2sigmapi)

    @property
    def n(self):
        return 3 * self.m + 2


class POAV2(POAV):
    """:class:`POAV` sub-class for POAV2 analysis."""

    def __init__(self, *args):
        super(POAV2, self).__init__(*args)

        bond_angles = self.bonds.angles
        bond_angle_pairs = self.bonds.bond_angle_pairs

        vi = []
        for bond, pair in zip(self.bonds, bond_angle_pairs):
            cosa = np.cos(bond_angles[np.in1d(self.bonds, pair, invert=True)])
            vi.append(cosa * bond.vector.unit_vector)

        self._v1 = vi[0]
        self._v2 = vi[1]
        self._v3 = vi[2]

        self.cosa12 = np.cos(bond_angles[0])
        self.cosa13 = np.cos(bond_angles[1])
        self.cosa23 = np.cos(bond_angles[2])

        self._T = -functools.reduce(operator.mul, np.cos(bond_angles), 1) * \
            self.T

    @property
    def n1(self):
        return -self.cosa23 / (self.cosa12 * self.cosa13)

    @property
    def n2(self):
        return -self.cosa13 / (self.cosa12 * self.cosa23)

    @property
    def n3(self):
        return -self.cosa12 / (self.cosa23 * self.cosa13)

    @property
    def m(self):
        s1 = 1 / (1 + self.n1)
        s2 = 1 / (1 + self.n2)
        s3 = 1 / (1 + self.n3)
        return 1 / sum([s1, s2, s3]) - 1


class POAVR(POAV):
    """:class:`POAV` sub-class for POAVR analysis."""

    def __init__(self, *args):
        super(POAVR, self).__init__(*args)

        vi = []
        for R, V in zip([self.R1, self.R2, self.R3],
                        [self.V1, self.V2, self.V3]):
            vi.append(R * V)

        self._v1 = vi[0]
        self._v2 = vi[1]
        self._v3 = vi[2]
        self._T = self.R1 * self.R2 * self.R3 * self.T
