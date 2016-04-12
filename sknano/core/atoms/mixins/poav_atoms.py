# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin classes for POAV analysis (:mod:`sknano.core.atoms.mixins.poav_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins.poav_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import functools
import operator
import warnings

from collections import OrderedDict

import numpy as np
np.seterr(all='warn')

from sknano.core import BaseClass, timethis
from sknano.core.math import Vectors, vector as vec

__all__ = ['POAV', 'POAV1', 'POAV2', 'POAVR',
           'POAVAtomMixin', 'POAVAtomsMixin']


class POAV(BaseClass):
    """Base class for POAV analysis.

    Parameters
    ----------
    atom : :class:`~sknano.core.atoms.Atom`
        :class:`~sknano.core.atoms.Atom` instance.

    """
    def __init__(self, atom):
        super().__init__()
        self.atom = atom
        self.angles = atom.angles
        self.bonds = atom.bonds

        self.bond1 = self.bonds[0].vector
        self.bond2 = self.bonds[1].vector
        self.bond3 = self.bonds[2].vector

        self.bond_angles = self.angles.angles
        self.angle_bond_pairs = self.angles.bond_pairs

        self.sigma_bond_angle12 = self.angles[0].angle
        self.sigma_bond_angle23 = self.angles[1].angle
        self.sigma_bond_angle31 = self.angles[2].angle

        self.cosa12 = np.cos(self.angles[0].angle)
        self.cosa23 = np.cos(self.angles[1].angle)
        self.cosa31 = np.cos(self.angles[2].angle)

        self._v1 = self.bond1
        self._v2 = self.bond2
        self._v3 = self.bond3

        self._pyramidalization_angles = None
        self._sigma_pi_angles = None
        self._misalignment_angles = None
        self.fmtstr = "{atom!r}"

    def __str__(self):
        fmtstr = '{}\n=====\n'.format(self.__class__.__name__)
        for k, v in self.POAVdict(rad2deg=True).items():
            fmtstr += '{}: {:.4f}\n'.format(k, v)
        return fmtstr

    @property
    def v1(self):
        """:class:`~sknano.core.math.Vector` :math:`\\mathbf{v}_1` \
            directed along the :math:`\\sigma`-orbital to the \
            nearest-neighbor :class:`~sknano.core.atoms.Atom` \
            in :class:`~sknano.core.atoms.mixins.Bond` 1."""
        return self._v1

    @property
    def v2(self):
        """:class:`~sknano.core.math.Vector` :math:`\\mathbf{v}_2` \
            directed along the :math:`\\sigma`-orbital to the \
            nearest-neighbor :class:`~sknano.core.atoms.Atom` \
            in :class:`~sknano.core.atoms.mixins.Bond` 2."""
        return self._v2

    @property
    def v3(self):
        """:class:`~sknano.core.math.Vector` :math:`\\mathbf{v}_3` \
            directed along the :math:`\\sigma`-orbital to the \
            nearest-neighbor :class:`~sknano.core.atoms.Atom` \
            in :class:`~sknano.core.atoms.mixins.Bond` 3."""
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
        """:class:`~sknano.core.atoms.mixins.Bond` 1 \
            :class:`~sknano.core.math.Vector` \
            :attr:`~sknano.core.math.Vector.length`.
        """
        return self.bond1.length

    @property
    def R2(self):
        """:class:`~sknano.core.atoms.mixins.Bond` 2 \
            :class:`~sknano.core.math.Vector` \
            :attr:`~sknano.core.math.Vector.length`.
        """
        return self.bond2.length

    @property
    def R3(self):
        """:class:`~sknano.core.atoms.mixins.Bond` 3 \
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

    def POAVdict(self, rad2deg=False):
        """Return dictionary of `POAV` class attributes."""
        bond_angles = self.bond_angles
        sigma_pi_angles = self.sigma_pi_angles
        pyramidalization_angles = self.pyramidalization_angles
        misalignment_angles = self.misalignment_angles

        if rad2deg:
            bond_angles = np.degrees(bond_angles)
            sigma_pi_angles = np.degrees(sigma_pi_angles)
            pyramidalization_angles = np.degrees(pyramidalization_angles)
            misalignment_angles = np.degrees(misalignment_angles)

        od = OrderedDict(
            [('bond1', self.bond1.length),
             ('bond2', self.bond2.length),
             ('bond3', self.bond3.length),
             ('sigma_bond_angle12', bond_angles[0]),
             ('sigma_bond_angle23', bond_angles[1]),
             ('sigma_bond_angle31', bond_angles[2]),
             ('sigma_pi_angle1', sigma_pi_angles[0]),
             ('sigma_pi_angle2', sigma_pi_angles[1]),
             ('sigma_pi_angle3', sigma_pi_angles[2]),
             ('pyramidalization_angle1', pyramidalization_angles[0]),
             ('pyramidalization_angle2', pyramidalization_angles[1]),
             ('pyramidalization_angle3', pyramidalization_angles[2]),
             ('misalignment_angle1', misalignment_angles[0]),
             ('misalignment_angle2', misalignment_angles[1]),
             ('misalignment_angle3', misalignment_angles[2]),
             ('T', self.T), ('H', self.H), ('A', self.A)])
        return od

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(atom=self.atom)


class POAV1(POAV):
    """:class:`POAV` sub-class for POAV1 analysis."""

    def __init__(self, *args):
        super().__init__(*args)

        self._v1 = self.V1
        self._v2 = self.V2
        self._v3 = self.V3

    @property
    def m(self):
        """:math:`s` character content of the :math:`\\pi`-orbital \
            (:math:`s^mp`) for :math:`sp^3` normalized hybridization."""
        cos2sigmapi = np.cos(np.mean(self.sigma_pi_angles)) ** 2
        return 2 * cos2sigmapi / (1 - 3 * cos2sigmapi)

    @property
    def n(self):
        """:math:`p` character content of the :math:`\\sigma`-orbitals \
            (:math:`sp^n`) for :math:`sp^3` normalized hybridization."""
        return 3 * self.m + 2

    def POAVdict(self, rad2deg=False):
        """Return dictionary of `POAV1` class attributes."""
        super_dict = super().POAVdict(rad2deg=rad2deg)
        super_dict.update([('m', self.m), ('n', self.n)])
        return super_dict


class POAV2(POAV):
    """:class:`POAV` sub-class for POAV2 analysis."""

    def __init__(self, *args):
        super().__init__(*args)
        vi = []
        for bond, pair in zip(self.bonds, self.angle_bond_pairs):
            # print('bond: {}'.format(bond))
            # print('pair: {}'.format(pair))
            # print(np.in1d(self.bonds, pair, invert=True))
            cosa = \
                np.cos(self.bond_angles[
                    np.in1d(self.bonds, pair, invert=True)])[0]
            # print('cosa: {}'.format(cosa))
            vi.append(cosa * bond.vector.unit_vector)
        # import sys
        # sys.exit(1)

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
                                 np.cos(self.bond_angles), 1) * super().T

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

    def POAVdict(self, rad2deg=False):
        """Return dictionary of `POAV2` class attributes."""
        super_dict = super().POAVdict(rad2deg=rad2deg)
        super_dict.update(
            [('m', self.m), ('n1', self.n1), ('n2', self.n2), ('n3', self.n3)])
        return super_dict


class POAVR(POAV):
    """:class:`POAV` sub-class for POAVR analysis."""

    def __init__(self, *args):
        super().__init__(*args)

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
        return self.R1 * self.R2 * self.R3 * super().T


class POAVAtomMixin:
    """Mixin class for :class:`POAV` analysis."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._POAV1 = None
        self._POAV2 = None
        self._POAVR = None

    @property
    def POAV1(self):
        """:class:`~sknano.core.atoms.mixins.POAV1` instance."""
        try:
            return self._POAV1
        except AttributeError:
            return None

    @POAV1.setter
    def POAV1(self, value):
        """Set :class:`~sknano.core.atoms.mixins.POAV1` instance."""
        if not isinstance(value, POAV1):
            raise TypeError('Expected a `POAV1` instance.')
        self._POAV1 = value

    @property
    def POAV2(self):
        """:class:`~sknano.core.atoms.mixins.POAV2` instance."""
        try:
            return self._POAV2
        except AttributeError:
            return None

    @POAV2.setter
    def POAV2(self, value):
        """Set :class:`~sknano.core.atoms.mixins.POAV2` instance."""
        if not isinstance(value, POAV2):
            raise TypeError('Expected a `POAV2` instance.')
        self._POAV2 = value

    @property
    def POAVR(self):
        """:class:`~sknano.core.atoms.mixins.POAVR` instance."""
        try:
            return self._POAVR
        except AttributeError:
            return None

    @POAVR.setter
    def POAVR(self, value):
        """Set :class:`~sknano.core.atoms.mixins.POAVR` instance."""
        if not isinstance(value, POAVR):
            raise TypeError('Expected a `POAVR` instance.')
        self._POAVR = value

    # def reset_attrs(self, poavs=False, **kwargs):
    #     """Reset :class:`POAVAtomMixin` attributes."""
    #     if poavs:
    #         [setattr(self, POAV, None) for POAV in
    #          ('_POAV1', '_POAV2', '_POAVR')]
    #     super().reset_attrs(**kwargs)


class POAVAtomsMixin:
    """Mixin class for POAV analysis."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.poavs = OrderedDict()

    @property
    def POAV1(self):
        """List of :class:`~sknano.core.atoms.mixins.POAVAtom` :class:`POAV1` \
            :attr:`~sknano.core.atoms.mixins.POAVAtom.POAV1` attribute."""
        return [atom.POAV1 for atom in self if atom.POAV1 is not None]

    @property
    def POAV2(self):
        """List of :class:`~sknano.core.atoms.mixins.POAVAtom` :class:`POAV2` \
            :attr:`~sknano.core.atoms.mixins.POAVAtom.POAV2` attribute."""
        return [atom.POAV2 for atom in self if atom.POAV2 is not None]

    @property
    def POAVR(self):
        """List of :class:`~sknano.core.atoms.mixins.POAVAtom` :class:`POAVR` \
            :attr:`~sknano.core.atoms.mixins.POAVAtom.POAVR` attribute."""
        return [atom.POAVR for atom in self if atom.POAVR is not None]

    @timethis
    def analyze_POAVs(self, **kwargs):
        """Compute `POAV1`, `POAV2`, `POAVR`."""
        if not self.neighbors_analyzed:
            self.update_neighbors(**kwargs)

        POAV_classes = {'POAV1': POAV1, 'POAV2': POAV2, 'POAVR': POAVR}

        for atom in self:
            # the central atom must have 3 bonds for POAV analysis.
            if atom.bonds.Nbonds == 3:
                for POAV_name, POAV_class in list(POAV_classes.items()):
                    setattr(atom, POAV_name, POAV_class(atom))

        for atom in self:
            # the central atom must have 3 bonds for POAV analysis.
            if atom.bonds.Nbonds == 3:
                for POAV_name in ('POAV1', 'POAV2', 'POAVR'):
                    POAV = getattr(atom, POAV_name)
                    sigma_pi_angles = []
                    pyramidalization_angles = []
                    misalignment_angles = []
                    for bond, NN in zip(atom.bonds, atom.NN):
                        # first compute the pyramidalization angle
                        sigma_pi_angle = vec.angle(POAV.Vpi, bond.vector)
                        if sigma_pi_angle < np.pi / 2:
                            sigma_pi_angle = np.pi - sigma_pi_angle
                        sigma_pi_angles.append(sigma_pi_angle)
                        pyramidalization_angles.append(
                            sigma_pi_angle - np.pi / 2)

                        # the bonded atom must have a POAV to compute the
                        # misalignment angles
                        if getattr(NN, POAV_name) is not None:
                            NN_POAV = getattr(NN, POAV_name)
                            # compute vector that is orthogonal to the plane
                            # defined by the bond vector and the POAV of the
                            # center atom.
                            nvec = vec.cross(bond.vector, POAV.Vpi)

                            # the misalignment angle is the angle between the
                            # nearest neighbor's POAV and the plane defined by
                            # the bond vector and the POAV of the center atom,
                            # which is pi/2 minus the angle between
                            # the NN POAV and the normal vector to the plane
                            # computed above.
                            misalignment_angles.append(np.abs(
                                np.pi / 2 - vec.angle(NN_POAV.Vpi, nvec)))
                        else:
                            misalignment_angles.append(np.nan)

                    POAV.pyramidalization_angles = pyramidalization_angles
                    POAV.misalignment_angles = misalignment_angles
                    POAV.sigma_pi_angles = sigma_pi_angles

    def compute_POAVs(self, **kwargs):
        """Alias for :meth:`~POAVAtomsMixin.analyze_POAVs`."""
        self.analyze_POAVs(**kwargs)

    def get_POAV_attr(self, POAV_class, attr):
        """Return list of :class:`~sknano.core.atoms.mixins.POAVAtom` :class:`POAV1` \
            :class:`POAV2` or :class:`POAVR` attribute.

        Parameters
        ----------
        POAV_class : {'POAV1', 'POAV2', 'POAVR'}
        attr : :class:`~python:str`

        Returns
        -------
        :class:`~python:tuple`

        """
        attr_values = []
        atom_ids = []
        attr_list = attr.split('.')
        for atom in self:
            POAV = getattr(atom, POAV_class)
            if POAV is not None:
                val = POAV
                for attr in attr_list:
                    val = getattr(val, attr)
                attr_values.append(val)
                atom_ids.append(atom.id)

        if len(attr_list) == 1 and attr_list[0] in POAV_vectors:
            attr_values = Vectors(attr_values)
        return attr_values, atom_ids

        # return [getattr(getattr(atom, POAV_class), attr) for atom in self
        #         if getattr(atom, POAV_class) is not None]

    def reset_attrs(self, poavs=False, **kwargs):
        """Reset the :class:`POAVAtomsMixin` class attributes, then call \
            parent class `reset_attrs` method."""
        if poavs:
            self.reset_poav_atoms_attrs()
        super().reset_attrs(**kwargs)

    def reset_poav_atoms_attrs(self):
        """Reset the :class:`POAVAtomsMixin` class attributes."""
        self.poavs = OrderedDict()
        [setattr(self, POAV, None) for POAV in ('_POAV1', '_POAV2', '_POAVR')]

    def update_attrs(self, poavs=False, **kwargs):
        """Update :class:`POAVAtomsMixin` class attributes."""
        super().update_attrs(**kwargs)
        if poavs:
            self.analyze_POAVs(**kwargs)


POAV_vectors = ('v1', 'v2', 'v3', 'vpi', 'V1', 'V2', 'V3', 'Vpi',
                'reciprocal_v1', 'reciprocal_v2', 'reciprocal_v3')
