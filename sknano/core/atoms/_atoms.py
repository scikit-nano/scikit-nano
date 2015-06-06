# -*- coding: utf-8 -*-
"""
==============================================================================
Base class for structure data atoms (:mod:`sknano.core.atoms._atoms`)
==============================================================================

.. currentmodule:: sknano.core.atoms._atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
from builtins import str
__docformat__ = 'restructuredtext en'

from operator import attrgetter

import numpy as np

from sknano.core import UserList
from ._atom import Atom

__all__ = ['Atoms']


class Atoms(UserList):
    """Base class for collection of `Atom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instance objects or an
        existing `Atoms` instance object.

    """
    def __init__(self, atoms=None, **kwargs):
        if atoms is not None and isinstance(atoms, str):
            atoms = [atoms]
        if not isinstance(atoms, type(self)) and isinstance(atoms, list):
            atoms = atoms[:]
            for i, atom in enumerate(atoms):
                try:
                    atoms[i] = self.__atom_class__(**atom.todict())
                except AttributeError:
                    atoms[i] = self.__atom_class__(atom, **kwargs)
        super().__init__(initlist=atoms)
        self.fmtstr = "{atoms!r}"
        self.kwargs = {}

    @property
    def __atom_class__(self):
        return Atom

    def __str__(self):
        """Return a nice string representation of `Atoms`."""
        return '\n'.join([str(atom) for atom in self])

    def __repr__(self):
        """Return canonical string representation of `Atoms`."""
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    def sort(self, key=attrgetter('element', 'Z', 'mass'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def fmtstr(self):
        return self._fmtstr

    @fmtstr.setter
    def fmtstr(self, value):
        self._fmtstr = value

    @property
    def Natoms(self):
        """Number of atoms in `Atoms`."""
        return len(self)

    @property
    def M(self):
        """Total mass of `Atoms`."""
        # return math.fsum(self.masses)
        return self.masses.sum()

    @property
    def elements(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`Atom.element`\ s \
            in `Atoms`."""
        return np.asarray([atom.element for atom in self])

    @property
    def masses(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`Atom.mass`\ s \
            in `Atoms`."""
        return np.asarray([atom.mass for atom in self])

    @property
    def symbols(self):
        """:class:`~numpy:numpy.ndarray` of :attr:`Atom.symbol`\ s \
            in `Atoms`."""
        return np.asarray([atom.symbol for atom in self])

    def filter(self, condition, invert=False):
        """Filter `Atoms` by `condition`.

        Parameters
        ----------
        condition : array_like, bool
            Boolean index array having same shape as the initial dimensions
            of the list of `Atoms` being indexed.
        invert : bool, optional
            If `True`, the boolean array `condition` is inverted element-wise.

        Returns
        -------
        filtered_atoms : `Atoms`
            If `invert` is `False`, return the elements where `condition`
            is `True`.

            If `invert` is `True`, return the elements where `~condition`
            (i.e., numpy.invert(condition)) is `True`.

        Examples
        --------
        An example using the structure data of a 10 nm `(10, 0)`
        `SWCNT`:

        >>> from sknano.generators import SWNTGenerator
        >>> swnt = SWNTGenerator(10, 0, Lz=10, fix_Lz=True).atoms
        >>> # select 'left', 'middle', 'right' atoms
        >>> latoms = swnt.filter(swnt.z <= 25)
        >>> matoms = swnt.filter((swnt.z < 75) & (swnt.z > 25))
        >>> ratoms = swnt.filter(swnt.z >= 75)
        >>> from pprint import pprint
        >>> pprint([getattr(atoms, 'bounds') for atoms in
        ...         (latoms, matoms, ratoms)])
        [Cuboid(pmin=Point([-3.914435, -3.914435, 0.0]), \
                pmax=Point([3.914435, 3.914435, 24.85])),
         Cuboid(pmin=Point([-3.914435, -3.914435, 25.56]), \
                pmax=Point([3.914435, 3.914435, 74.55])),
         Cuboid(pmin=Point([-3.914435, -3.914435, 75.97]), \
                pmax=Point([3.914435, 3.914435, 100.11]))]
        >>> latoms.Natoms + matoms.Natoms + ratoms.Natoms == swnt.Natoms
        True

        """
        if invert:
            condition = ~condition
        return self.__class__(atoms=np.asarray(self)[condition].tolist(),
                              **self.kwargs)

    def get_atoms(self, asarray=False):
        """Return list of `Atoms`.

        Parameters
        ----------
        asarray : bool, optional

        Returns
        -------
        sequence or ndarray

        """
        if asarray:
            return np.asarray(self.data)
        else:
            return self.data

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

        """
        try:
            return np.asarray([getattr(atom, attr) for atom in self])
        except AttributeError:
            return None

    def mapatomattr(self, attr, from_attr, attrmap):
        """Set/update atom attribute from another atom attribute with dict.

        Parameters
        ----------
        attr, from_attr : :class:`python:str`
        attrmap : :class:`python:dict`

        Examples
        --------
        Suppose you have an `XAtoms` instance named ``atoms`` that has
        `XAtom` instances of two atom types `1` and `2` and we want to set
        all `XAtom`\ s with `type=1` to Nitrogen and all `XAtom`\ s with
        `type=2` to Argon.

        We'd call this method like so::

        >>> atoms.mapatomattr('element', 'type', {1: 'N', 2: 'Ar'})

        """
        [setattr(atom, attr, attrmap[getattr(atom, from_attr)])
         for atom in self if getattr(atom, from_attr) is not None]

    def rotate(self, **kwargs):
        assert not hasattr(super(), 'rotate')

    def translate(self, *args, **kwargs):
        assert not hasattr(super(), 'translate')

    def todict(self):
        return dict(atoms=self.data)
