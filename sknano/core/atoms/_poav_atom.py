# -*- coding: utf-8 -*-
"""
===============================================================================
Atom sub-class for POAV analysis (:mod:`sknano.core.atoms._poav_atom`)
===============================================================================

An `Atom` sub-class for POAV analysis.

.. currentmodule:: sknano.core.atoms._poav_atom

"""
from __future__ import absolute_import, division, print_function
from six.moves import zip
__docformat__ = 'restructuredtext en'

import functools
import operator
import warnings

import numpy as np
np.seterr(all='warn')

from sknano.core.math import vector as vec

from ._kdtree_atom import KDTAtom

__all__ = ['POAV', 'POAV1', 'POAV2', 'POAVR', 'POAVAtomMixin', 'POAVAtom']


class POAV(object):
    """Base class for POAV analysis.

    Parameters
    ----------
    sigma_bonds : `~sknano.core.atoms.Bonds`
        `~sknano.core.atoms.Bonds` instance.

    """
    def __init__(self, sigma_bonds):
        self.bonds = sigma_bonds
        self.bond1 = self.bonds[0].vector
        self.bond2 = self.bonds[1].vector
        self.bond3 = self.bonds[2].vector

        self.bond_angles = self.bonds.angles
        self.bond_angle_pairs = self.bonds.bond_angle_pairs

        self.sigma_bond_angle12 = self.bond_angles[0]
        self.sigma_bond_angle23 = self.bond_angles[1]
        self.sigma_bond_angle31 = self.bond_angles[2]

        self.cosa12 = np.cos(self.bond_angles[0])
        self.cosa23 = np.cos(self.bond_angles[1])
        self.cosa31 = np.cos(self.bond_angles[2])

        self._v1 = self.bond1
        self._v2 = self.bond2
        self._v3 = self.bond3

        self._pyramidalization_angles = None
        self._sigma_pi_angles = None
        self._misalignment_angles = None

    @property
    def v1(self):
        """:class:`~sknano.core.math.Vector` :math:`\\mathbf{v}_1` \
            directed along the :math:`\\sigma`-orbital to the \
            nearest-neighbor :class:`~sknano.core.atoms.Atom` \
            in :class:`~sknano.core.atoms.Bond` 1."""
        return self._v1

    @property
    def v2(self):
        """:class:`~sknano.core.math.Vector` :math:`\\mathbf{v}_2` \
            directed along the :math:`\\sigma`-orbital to the \
            nearest-neighbor :class:`~sknano.core.atoms.Atom` \
            in :class:`~sknano.core.atoms.Bond` 2."""
        return self._v2

    @property
    def v3(self):
        """:class:`~sknano.core.math.Vector` :math:`\\mathbf{v}_3` \
            directed along the :math:`\\sigma`-orbital to the \
            nearest-neighbor :class:`~sknano.core.atoms.Atom` \
            in :class:`~sknano.core.atoms.Bond` 3."""
        return self._v3

    @property
    def Vv1v2v3(self):
        """Volume of the parallelepiped defined by \
            :class:`~sknano.core.math.Vector`\ s `v1`, `v2`, and `v3`.

        Computes the scalar triple product of vectors :math:`\\mathbf{v}_1`,
        :math:`\\mathbf{v}_2`, and :math:`\\mathbf{v}_3`:

        .. math::

           V_{v_1v_2v_3} =
           |\\mathbf{v}_1\\cdot(\\mathbf{v}_2\\times\\mathbf{v}_3)|

        """
        return np.abs(vec.scalar_triple_product(self.v1, self.v2, self.v3))

    @property
    def vpi(self):
        """General :math:`\\pi`-orbital axis vector \
            (:math:`\\mathbf{v}_{\\pi}`) formed by the \
            terminii of :class:`~sknano.core.math.Vector`\ s \
            :class:`~sknano.core.math.Vector`\ s `v1`, `v2`, and `v3`.

        .. math::

           \\mathbf{v}_{\\pi} =
           \\mathbf{v}_1 + \\mathbf{v}_2\\ + \\mathbf{v}_3

        """
        return self.reciprocal_v1 + self.reciprocal_v2 + self.reciprocal_v3

    @property
    def Vpi(self):
        """:math:`\\mathbf{v}_{\\pi}` unit :class:`~sknano.core.math.Vector`

        Returns the :math:`\\pi`-orbital axis vector
        (:math:`\\mathbf{v}_{\\pi}`) unit vector.

        .. math::

           \\mathbf{V}_{\\pi} =
           \\frac{\\mathbf{v}_{\\pi}}{|\\mathbf{v}_{\\pi}|}

        """
        return self.vpi.unit_vector

    @property
    def reciprocal_v1(self):
        """Reciprocal :class:`~sknano.core.math.Vector` \
            :math:`\\mathbf{v}_1^{*}`.

        Defined as:

        .. math::

           \\mathbf{v}_1^{*} =
           \\frac{\\mathbf{v}_2\\times\\mathbf{v}_3}
           {|\\mathbf{v}_1\\cdot(\\mathbf{v}_2\\times\\mathbf{v}_3)|}

        """
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                return vec.cross(self.v2, self.v3) / self.Vv1v2v3
            except Warning:
                return vec.cross(self.v2, self.v3)

    @property
    def reciprocal_v2(self):
        """Reciprocal :class:`~sknano.core.math.Vector` \
            :math:`\\mathbf{v}_2^{*}`.

        Defined as:

        .. math::

           \\mathbf{v}_2^{*} =
           \\frac{\\mathbf{v}_3\\times\\mathbf{v}_1}
           {|\\mathbf{v}_1\\cdot(\\mathbf{v}_2\\times\\mathbf{v}_3)|}

        """
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                return vec.cross(self.v3, self.v1) / self.Vv1v2v3
            except Warning:
                return vec.cross(self.v3, self.v1)

    @property
    def reciprocal_v3(self):
        """Reciprocal :class:`~sknano.core.math.Vector` \
            :math:`\\mathbf{v}_3^{*}`.

        Defined as:

        .. math::

           \\mathbf{v}_3^{*} =
           \\frac{\\mathbf{v}_1\\times\\mathbf{v}_2}
           {|\\mathbf{v}_1\\cdot(\\mathbf{v}_2\\times\\mathbf{v}_3)|}

        """
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                return vec.cross(self.v1, self.v2) / self.Vv1v2v3
            except Warning:
                return vec.cross(self.v1, self.v2)

    @property
    def V1(self):
        """:math:`\\mathbf{v}_1` unit :class:`~sknano.core.math.Vector`

        .. math::
           \\mathbf{V}_1\\equiv\\frac{\\mathbf{v}_1}{|\\mathbf{v}_1|}

        """
        return self.bond1.unit_vector

    @property
    def V2(self):
        """:math:`\\mathbf{v}_2` unit :class:`~sknano.core.math.Vector`

        .. math::
           \\mathbf{V}_2\\equiv\\frac{\\mathbf{v}_2}{|\\mathbf{v}_2|}

        """
        return self.bond2.unit_vector

    @property
    def V3(self):
        """:math:`\\mathbf{v}_3` unit :class:`~sknano.core.math.Vector`

        .. math::
           \\mathbf{V}_3\\equiv\\frac{\\mathbf{v}_3}{|\\mathbf{v}_3|}

        """
        return self.bond3.unit_vector

    @property
    def R1(self):
        """:class:`~sknano.core.atoms.Bond` 1 \
            :class:`~sknano.core.math.Vector` \
            :attr:`~sknano.core.math.Vector.length`.
        """
        return self.bond1.length

    @property
    def R2(self):
        """:class:`~sknano.core.atoms.Bond` 2 \
            :class:`~sknano.core.math.Vector` \
            :attr:`~sknano.core.math.Vector.length`.
        """
        return self.bond2.length

    @property
    def R3(self):
        """:class:`~sknano.core.atoms.Bond` 3 \
            :class:`~sknano.core.math.Vector` \
            :attr:`~sknano.core.math.Vector.length`.
        """
        return self.bond3.length

    @property
    def t(self):
        """:math:`\\frac{1}{6}` the volume of the tetrahedron defined by \
            :class:`~sknano.core.math.Vector`\ s `v1`, `v2`, and `v3`.

        .. math::
           t =
           \\frac{|\\mathbf{v}_1\\cdot(\\mathbf{v}_2\\times\\mathbf{v}_3)|}{6}

        """
        return self.Vv1v2v3 / 6

    @property
    def T(self):
        """:math:`\\frac{1}{6}` the volume of the tetrahedron defined by \
            :class:`~sknano.core.math.Vector`\ s `V1`, `V2`, and `V3`.

        .. math::
           T =
           \\frac{|\\mathbf{V}_1\\cdot(\\mathbf{V}_2\\times\\mathbf{V}_3)|}{6}

        """
        return np.abs(vec.scalar_triple_product(self.V1, self.V2, self.V3) / 6)

    @property
    def A(self):
        """Magnitude of :math:`\\mathbf{v}_{\\pi}`."""
        return self.vpi.magnitude

    @property
    def H(self):
        """Altitude of tetrahedron."""
        return 3 * self.T / self.A

    @property
    def sigma_pi_angles(self):
        """List of :math:`\\theta_{\\sigma-\\pi}` angles."""
        return self._sigma_pi_angles

    @sigma_pi_angles.setter
    def sigma_pi_angles(self, value):
        """Set list of :math:`\\theta_{\\sigma-\\pi}` angles."""
        if not isinstance(value, list):
            raise TypeError('Expected a list')
        self._sigma_pi_angles = value

    @property
    def pyramidalization_angles(self):
        """List of pyramidalization :math:`\\theta_{P}` angles."""
        return self._pyramidalization_angles

    @pyramidalization_angles.setter
    def pyramidalization_angles(self, value):
        """Set list of :math:`\\theta_{P}` angles."""
        if not isinstance(value, list):
            raise TypeError('Expected a list')
        self._pyramidalization_angles = value

    @property
    def misalignment_angles(self):
        """List of misalignment :math:`\\phi_{i}` angles."""
        return self._misalignment_angles

    @misalignment_angles.setter
    def misalignment_angles(self, value):
        """Set list of :math:`\\phi` angles."""
        if not isinstance(value, list):
            raise TypeError('Expected a list')
        self._misalignment_angles = value

    def todict(self, rad2deg=False):
        """Return dictionary of `POAV` class attributes."""
        sigma_pi_angles = self.sigma_pi_angles
        pyramidalization_angles = self.pyramidalization_angles
        misalignment_angles = self.misalignment_angles

        if rad2deg:
            sigma_pi_angles = np.degrees(sigma_pi_angles)
            pyramidalization_angles = np.degrees(pyramidalization_angles)
            misalignment_angles = np.degrees(misalignment_angles)

        return dict(bond1=self.bond1.length,
                    bond2=self.bond2.length,
                    bond3=self.bond3.length,
                    sigma_bond_angle12=self.sigma_bond_angle12,
                    sigma_bond_angle23=self.sigma_bond_angle23,
                    sigma_bond_angle31=self.sigma_bond_angle31,
                    sigma_pi_angle1=sigma_pi_angles[0],
                    sigma_pi_angle2=sigma_pi_angles[1],
                    sigma_pi_angle3=sigma_pi_angles[2],
                    pyramidalization_angle1=pyramidalization_angles[0],
                    pyramidalization_angle2=pyramidalization_angles[1],
                    pyramidalization_angle3=pyramidalization_angles[2],
                    misalignment_angle1=misalignment_angles[0],
                    misalignment_angle2=misalignment_angles[1],
                    misalignment_angle3=misalignment_angles[2],
                    T=self.T, H=self.H, A=self.A)


class POAV1(POAV):
    """:class:`POAV` sub-class for POAV1 analysis."""

    def __init__(self, *args):
        super(POAV1, self).__init__(*args)

        self._v1 = self.V1
        self._v2 = self.V2
        self._v3 = self.V3

    @property
    def m(self):
        """:math:`s` character content of the :math:`\\pi`-orbital \
            (:math:`s^mp`) for :math:`sp^3` normalized hybridization."""
        cos2sigmapi = np.cos(np.mean(self.sigma_pi_angles))**2
        return 2 * cos2sigmapi / (1 - 3 * cos2sigmapi)

    @property
    def n(self):
        """:math:`p` character content of the :math:`\\sigma`-orbitals \
            (:math:`sp^n`) for :math:`sp^3` normalized hybridization."""
        return 3 * self.m + 2

    def todict(self, rad2deg=False):
        """Return dictionary of `POAV1` class attributes."""
        super_dict = super(POAV1, self).todict(rad2deg=rad2deg)
        super_dict.update(dict(m=self.m, n=self.n))
        return super_dict


class POAV2(POAV):
    """:class:`POAV` sub-class for POAV2 analysis."""

    def __init__(self, *args):
        super(POAV2, self).__init__(*args)

        vi = []
        for bond, pair in zip(self.bonds, self.bond_angle_pairs):
            cosa = \
                np.cos(self.bond_angles[
                    np.in1d(self.bonds, pair, invert=True)])
            vi.append(cosa * bond.vector.unit_vector)

        self._v1 = vi[0]
        self._v2 = vi[1]
        self._v3 = vi[2]

    @property
    def T(self):
        """:math:`\\frac{1}{6}` the volume of the tetrahedron defined by \
            :class:`~sknano.core.math.Vector`\ s `V1`, `V2`, and `V3`.

        .. math::
           T =
           \\cos\\theta_{12}\\cos\\theta_{23}\\cos\\theta_{31}\\times
           \\frac{|\\mathbf{V}_1\\cdot(\\mathbf{V}_2\\times\\mathbf{V}_3)|}{6}

        """
        return -functools.reduce(operator.mul,
                                 np.cos(self.bonds.angles), 1) * \
            super(POAV2, self).T

    @property
    def n1(self):
        """:math:`p` character content of the :math:`\\sigma`-orbital \
            hybridization for :math:`\\sigma_1` bond."""
        return -self.cosa23 / (self.cosa12 * self.cosa31)

    @property
    def n2(self):
        """:math:`p` character content of the :math:`\\sigma`-orbital \
            hybridization for :math:`\\sigma_2` bond."""
        return -self.cosa31 / (self.cosa12 * self.cosa23)

    @property
    def n3(self):
        """:math:`p` character content of the :math:`\\sigma`-orbital \
            hybridization for :math:`\\sigma_3` bond."""
        return -self.cosa12 / (self.cosa31 * self.cosa23)

    @property
    def m(self):
        """:math:`s` character content of the :math:`\\pi`-orbital \
            (:math:`s^mp`) for :math:`sp^3` normalized hybridization."""
        s1 = 1 / (1 + self.n1)
        s2 = 1 / (1 + self.n2)
        s3 = 1 / (1 + self.n3)
        return 1 / sum([s1, s2, s3]) - 1

    def todict(self, rad2deg=False):
        """Return dictionary of `POAV2` class attributes."""
        super_dict = super(POAV2, self).todict(rad2deg=rad2deg)
        super_dict.update(dict(m=self.m, n1=self.n1, n2=self.n2, n3=self.n3))
        return super_dict


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

    @property
    def T(self):
        """:math:`\\frac{1}{6}` the volume of the tetrahedron defined by \
            :class:`~sknano.core.math.Vector`\ s `V1`, `V2`, and `V3`.

        .. math::
           T =
           R_1 R_2 R_3 \\times
           \\frac{|\\mathbf{V}_1\\cdot(\\mathbf{V}_2\\times\\mathbf{V}_3)|}{6}

        """
        return self.R1 * self.R2 * self.R3 * super(POAVR, self).T


class POAVAtomMixin(object):
    """Mixin class for `POAVAtom`."""
    @property
    def POAV1(self):
        """:class:`~sknano.utils.analysis.POAV1` instance."""
        return self._POAV1

    @POAV1.setter
    def POAV1(self, value):
        """Set :class:`~sknano.utils.analysis.POAV1` instance."""
        if not isinstance(value, POAV1):
            raise TypeError('Expected a `POAV1` instance.')
        self._POAV1 = value

    @property
    def POAV2(self):
        """:class:`~sknano.utils.analysis.POAV2` instance."""
        return self._POAV2

    @POAV2.setter
    def POAV2(self, value):
        """Set :class:`~sknano.utils.analysis.POAV2` instance."""
        if not isinstance(value, POAV2):
            raise TypeError('Expected a `POAV2` instance.')
        self._POAV2 = value

    @property
    def POAVR(self):
        """:class:`~sknano.utils.analysis.POAVR` instance."""
        return self._POAVR

    @POAVR.setter
    def POAVR(self, value):
        """Set :class:`~sknano.utils.analysis.POAVR` instance."""
        if not isinstance(value, POAVR):
            raise TypeError('Expected a `POAVR` instance.')
        self._POAVR = value


class POAVAtom(POAVAtomMixin, KDTAtom):
    """An `Atom` sub-class for POAV analysis."""
    _atomattrs = KDTAtom._atomattrs + ['POAV1', 'POAV2', 'POAVR']

    def __init__(self, **kwargs):
        super(POAVAtom, self).__init__(**kwargs)
        self._POAV1 = None
        self._POAV2 = None
        self._POAVR = None
