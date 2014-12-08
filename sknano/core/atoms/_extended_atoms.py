# -*- coding: utf-8 -*-
"""
===============================================================================
Extended Atoms class feature set (:mod:`sknano.core.atoms._extended_atoms`)
===============================================================================

An "eXtended" `Atoms` class for structure analysis.

.. currentmodule:: sknano.core.atoms._extended_atoms

"""
from __future__ import absolute_import, division, print_function
from six.moves import zip
__docformat__ = 'restructuredtext en'

from collections import OrderedDict
from operator import attrgetter

import numpy as np

from sknano.core import xyz
from sknano.core.math import Vector, transformation_matrix
from sknano.utils.geometric_shapes import Cuboid  # , Rectangle
from ._atoms import Atoms

__all__ = ['XAtoms']


class XAtoms(Atoms):
    """An eXtended `Atoms` class.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.XAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `XAtoms`}, optional
        if not `None`, then a list of `XAtom` instance objects or an
        existing `XAtoms` instance object.
    copylist : bool, optional
        perform shallow copy of atoms list
    deepcopy : bool, optional
        perform deepcopy of atoms list

    """
    _atomattrs = Atoms._atomattrs + \
        ['atomID', 'moleculeID', 'atomtype', 'q', 'dr',
         'r', 'x', 'y', 'z', 'v', 'vx', 'vy', 'vz',
         'f', 'fx', 'fy', 'fz', 'n', 'nx', 'ny', 'nz']

    def __init__(self, atoms=None, copylist=True, deepcopy=False):

        super(XAtoms, self).__init__(atoms=atoms,
                                     copylist=copylist,
                                     deepcopy=deepcopy)
        self._atomtypes = {}

    def sort(self, key=None, reverse=False):
        if key is None:
            self.data.sort(key=attrgetter('element', 'Z', 'atomtype',
                                          'moleculeID', 'atomID', 'z'),
                           reverse=reverse)
        else:
            self.data.sort(key=key, reverse=reverse)

    @property
    def CM(self):
        """Center-of-Mass coordinates of `Atoms`.

        Computes the position vector of the center-of-mass coordinates:

        .. math::

           \\mathbf{R}_{CM} = \\frac{1}{M}\\sum_{i=1}^{N_{\\mathrm{atoms}}}
           m_i\\mathbf{r}_i

        Returns
        -------
        CM : :class:`~sknano.core.math.Vector`
            The position vector of the center of mass coordinates.

        """
        masses = np.asarray([self.masses])
        coords = self.coords
        MxR = masses.T * coords
        CM = Vector(np.sum(MxR, axis=0) / np.sum(masses))
        CM.rezero()
        return CM

    @property
    def Ntypes(self):
        """Number of :attr:`~XAtoms.atomtypes`."""
        return len(list(self.atomtypes.keys()))

    @property
    def centroid(self):
        """Centroid of `Atoms`.

        Computes the position vector of the centroid of the `Atoms`
        coordinates.

        .. math::
           \\mathbf{C} =
           \\frac{\\sum_{i=1}^{N_{\\mathrm{atoms}}}
           m_i\\mathbf{r}_i}{\\sum_{i=1}^{N_{\\mathrm{atoms}}}m_i}

        Returns
        -------
        C : :class:`~sknano.core.math.Vector`
            The position vector of the centroid coordinates.
        """
        C = Vector(np.mean(self.coords, axis=0))
        C.rezero()
        return C

    @property
    def bounds(self):
        """Bounds of `Atoms`.

        Returns
        -------
        :class:`~sknano.utils.geometric_shapes.Cuboid`"""
        return Cuboid(pmin=[self.x.min(), self.y.min(), self.z.min()],
                      pmax=[self.x.max(), self.y.max(), self.z.max()])

    @property
    def coords(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s \
            :math:`x, y, z` coordinates."""
        return np.asarray([atom.r for atom in self])

    @property
    def x(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`x` coordinates."""
        return self.coords[:,0]

    @property
    def y(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`y` coordinates."""
        return self.coords[:,1]

    @property
    def z(self):
        """:class:`~numpy:numpy.ndarray` of `Atom`\ s :math:`z` coordinates."""
        return self.coords[:,2]

    @property
    def inertia_tensor(self):
        """Return the inertia tensor."""
        Ixx = (self.masses * (self.y**2 + self.z**2)).sum()
        Iyy = (self.masses * (self.x**2 + self.z**2)).sum()
        Izz = (self.masses * (self.x**2 + self.y**2)).sum()
        Ixy = Iyx = (-self.masses * self.x * self.y).sum()
        Ixz = Izx = (-self.masses * self.x * self.z).sum()
        Iyz = Izy = (-self.masses * self.y * self.z).sum()
        return np.array([[Ixx, Ixy, Ixz], [Iyx, Iyy, Iyz], [Izx, Izy, Izz]])

    @property
    def atomtypes(self):
        """:attr:`~XAtoms.atomtypes` :class:`python:dict`."""
        self._update_atomtypes()
        return self._atomtypes

    def _update_atomtypes(self):
        for atom in self:
            if atom.atomtype not in self._atomtypes:
                self._atomtypes[atom.atomtype] = {}
                self._atomtypes[atom.atomtype]['mass'] = atom.m
                self._atomtypes[atom.atomtype]['q'] = atom.q

    @property
    def atom_ids(self):
        """Return array of `XAtom` IDs."""
        if len(set([atom.atomID for atom in self])) != self.Natoms:
            self.assign_unique_ids()
        return np.asarray([atom.atomID for atom in self])

    @property
    def charges(self):
        """Return array of `XAtom` charges."""
        return np.asarray([atom.q for atom in self])

    @property
    def q(self):
        """Return the total net charge of `XAtoms`."""
        return self.charges.sum()

    @property
    def velocities(self):
        """Return array of `XAtom` velocities."""
        return np.asarray([atom.v for atom in self])

    def add_atomtype(self, atom):
        """Add atom type to :attr:`~XAtoms.atomtypes`.

        Parameters
        ----------
        atom : :class:`~sknano.core.atoms.XAtom`
            A :class:`~sknano.core.atoms.XAtom` instance.

        """
        if atom.atomtype not in self._atomtypes:
            self._atomtypes[atom.atomtype] = {}
            self._atomtypes[atom.atomtype]['mass'] = atom.m
            self._atomtypes[atom.atomtype]['q'] = atom.q

    def add_atomtypes(self, atoms=None):
        """Add atomtype for each atom in atoms to atomtypes dictionary.
        Parameters
        ----------
        atoms : sequence
            a list of `XAtom` object instances

        """
        try:
            [self.add_atomtype(atom) for atom in atoms]
        except TypeError:
            print('Expected an iterable sequence of `XAtom` objects.')

    def assign_unique_ids(self, starting_id=1):
        """Assign unique ID to each `XAtom` in `XAtoms`."""
        [setattr(atom, 'atomID', i)
         for i, atom in enumerate(self, start=starting_id)]

    def center_CM(self, axes=None):
        """Center atoms on CM coordinates."""
        dr = -self.CM
        self.translate(dr)

    def clip_bounds(self, region, center_before_clipping=False):
        """Remove atoms outside the given limits along given dimension.

        Parameters
        ----------
        region : :class:`~sknano.utils.geometric_shapes.`GeometricRegion`

        """
        CM0 = None
        if center_before_clipping:
            CM0 = self.CM
            self.translate(-CM0)

        self.data = \
            np.asarray(self)[np.logical_and(
                np.logical_and(
                    self.x <= region.limits['x']['max'],
                    np.logical_and(
                        self.y <= region.limits['y']['max'],
                        self.z <= region.limits['z']['max'])),
                np.logical_and(
                    self.x >= region.limits['x']['min'],
                    np.logical_and(
                        self.y >= region.limits['y']['min'],
                        self.z >= region.limits['z']['min'])))].tolist()

        #for dim, limits in region.limits.iteritems():
        #    atoms = atoms[np.where(getattr(self, dim) <= limits['max'])]
        #    atoms = atoms[np.where(getattr(self, dim) >= limits['min'])]
        #    self = atoms.tolist()

        if CM0 is not None:
            self.translate(CM0)

    def filter_ids(self, atom_ids, invert=False):
        """Return `Atoms` by :attr:`XAtoms.atom_ids` in `atom_ids`.

        Parameters
        ----------
        filter_ids : array_like
        invert : bool, optional

        Returns
        -------
        filtered_atoms : `Atoms`
            An instance of `Atoms` (sub)class.

        """
        filtered_atoms = \
            np.asarray(self)[np.in1d(
                self.atom_ids, atom_ids, invert=invert).nonzero()].tolist()
        return self.__class__(atoms=filtered_atoms)

    def get_atom(self, atomID):
        """Get `XAtom` with :attr:`Xatom.atomID` == `atomID`.

        Parameters
        ----------
        atomID : int

        Returns
        -------
        `XAtom` instance with :attr:`XAtom.atomID` == `atomID` or
        `None` if no :attr:`XAtom.atomID` == `atomID`.

        """
        try:
            return self[np.where(self.atom_ids == atomID)[0]]
        except TypeError:
            print('No atom with atomID = {}'.format(atomID))
            return None

    def get_coords(self, asdict=False):
        """Return atom coords.

        Parameters
        ----------
        asdict : bool, optional

        Returns
        -------
        coords : :py:class:`python:~collections.OrderedDict` or ndarray

        """
        coords = self.coords
        if asdict:
            return OrderedDict(list(zip(xyz, coords.T)))
        else:
            return coords

    def getatomattr(self, attr):
        """Get :class:`~numpy:numpy.ndarray` of atom attributes `attr`.

        Parameters
        ----------
        attr : str
            Name of attribute to pass to `getattr` from each `XAtom` in
            `XAtoms`.

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        Raises
        ------
        `ValueError`
            if `attr` not in list of class atomattrs.

        """
        if attr not in self._atomattrs:
            errmsg = '{} not in list of allowed attributes:\n{}'
            raise ValueError(errmsg.format(attr, self._atomattrs))
        return np.asarray([getattr(atom, attr) for atom in self])

    def mapatomattr(self, attr, from_attr, attrmap):
        """Set/update atom attribute from another atom attribute with dict.

        Parameters
        ----------
        attr, from_attr : :class:`python:str`
        attrmap : :class:`python:dict`

        Examples
        --------
        Suppose you have an `XAtoms` instance named ``atoms`` that has
        `XAtom` instances of two atomtypes `1` and `2` and we want to set
        all `XAtom`\ s with `atomtype=1` to Nitrogen and all `XAtom`\ s with
        `atomtype=2` to Argon.

        We'd call this method like so::

        >>> atoms.mapatomattr('element', 'atomtype', {1: 'N', 2: 'Ar'})

        """
        [setattr(atom, attr, attrmap[getattr(atom, from_attr)])
         for atom in self if getattr(atom, from_attr) is not None]

    def rezero_coords(self, epsilon=1.0e-10):
        """Alias for :meth:`Atoms.rezero`."""
        self.rezero(epsilon=epsilon)

    def rezero_xyz(self, epsilon=1.0e-10):
        """Alias for :meth:`Atoms.rezero`."""
        self.rezero(epsilon=epsilon)

    def rezero(self, epsilon=1.0e-10):
        """Set really really small coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        [atom.rezero(epsilon=epsilon) for atom in self]

    def rotate(self, angle=None, rot_axis=None, anchor_point=None,
               rot_point=None, from_vector=None, to_vector=None,
               deg2rad=False, transform_matrix=None, verbose=False):
        """Rotate `Atom` position vectors.

        Parameters
        ----------
        angle : float
        rot_axis : :class:`~sknano.core.math.Vector`, optional
        anchor_point : :class:`~sknano.core.math.Point`, optional
        rot_point : :class:`~sknano.core.math.Point`, optional
        from_vector, to_vector : :class:`~sknano.core.math.Vector`, optional
        deg2rad : bool, optional
        transform_matrix : :class:`~numpy:numpy.ndarray`

        """
        if transform_matrix is None:
            transform_matrix = \
                transformation_matrix(angle=angle, rot_axis=rot_axis,
                                      anchor_point=anchor_point,
                                      rot_point=rot_point,
                                      from_vector=from_vector,
                                      to_vector=to_vector, deg2rad=deg2rad,
                                      verbose=verbose)
        [atom.rotate(transform_matrix=transform_matrix) for atom in self]

    def select(self, cmd):
        pass

    def select_within(self, volume):
        pass

    def translate(self, t, fix_anchor_points=True):
        """Translate `Atom` position vectors by :class:`Vector` `t`.

        Parameters
        ----------
        t : :class:`Vector`
        fix_anchor_points : bool, optional

        """
        [atom.translate(t, fix_anchor_point=fix_anchor_points)
         for atom in self]
