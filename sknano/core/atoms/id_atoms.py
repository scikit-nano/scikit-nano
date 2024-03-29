# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes with id attributes (:mod:`sknano.core.atoms.id_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms.id_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

from operator import attrgetter
import numbers

import numpy as np

from .atoms import Atom, Atoms

__all__ = ['IDAtom', 'IDAtoms']


class IDAtom(Atom):
    """An `Atom` sub-class with id attributes.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    id : int, optional
        atom ID
    mol : int, optional
        molecule ID
    """

    def __init__(self, *args, id=0, serial=0, mol=0, **kwargs):

        if 'atomID' in kwargs:
            id = kwargs['atomID']
            del kwargs['atomID']

        if 'moleculeID' in kwargs:
            mol = kwargs['moleculeID']
            del kwargs['moleculeID']

        super().__init__(*args, **kwargs)

        if serial != id and id == 0:
            id = serial
        self.id = id
        self.mol = mol
        self.fmtstr = super().fmtstr + ", id={id!r}, mol={mol!r}"

    @property
    def __atoms_class__(self):
        return IDAtoms

    # def __eq__(self, other):
    #     return self.id == other.id and self.mol == other.mol and \
    #         super().__eq__(other)

    # def __lt__(self, other):
    #     return ((self.id < other.id and self.mol <= other.mol and
    #              super().__le__(other))
    #             or (self.id <= other.id and self.mol < other.mol and
    #                 super().__le__(other))
    #             or (self.id <= other.id and self.mol <= other.mol and
    #                 super().__lt__(other)))

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.id == other.id and super().__eq__(other)

    def __le__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.id > other.id or not super().__le__(other):
            return False
        return True

    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.id >= other.id or not super().__lt__(other):
            return False
        return True

    def __ge__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.id < other.id or not super().__ge__(other):
            return False
        return True

    def __gt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.id <= other.id or not super().__gt__(other):
            return False
        return True

    def __dir__(self):
        attrs = super().__dir__()
        attrs.extend(['id', 'mol'])
        return attrs

    @property
    def id(self):
        """:attr:`~IDAtom.id`."""
        return self._id

    @id.setter
    def id(self, value):
        """Set atom id.

        Parameters
        ----------
        value : :class:`~python:int`
            atom ID

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._id = int(value)

    @property
    def mol(self):
        """:attr:`~IDAtom.mol`."""
        return self._mol

    @mol.setter
    def mol(self, value):
        """Set :attr:`~IDAtom.mol`.

        Parameters
        ----------
        value : int
            molecule ID

        """
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected a number')
        self._mol = int(value)

    @property
    def molid(self):
        """An alias for :attr:`~IDAtom.mol`."""
        return self.mol

    @molid.setter
    def molid(self, value):
        self.mol = value

    @property
    def atomID(self):
        """Alias for :attr:`~IDAtom.id`."""
        return self.id

    @atomID.setter
    def atomID(self, value):
        self.id = value

    @property
    def serial(self):
        """Alias for :attr:`~IDAtom.id`."""
        return self.id

    @serial.setter
    def serial(self, value):
        self.id = value

    @property
    def moleculeID(self):
        """Alias for :attr:`~IDAtom.mol`."""
        return self.mol

    @moleculeID.setter
    def moleculeID(self, value):
        self.mol = value

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        super_dict = super().todict()
        super_dict.update(dict(id=self.id, mol=self.mol))
        return super_dict


class IDAtoms(Atoms):
    """An `Atoms` sub-class for `IDAtom`\ s.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.IDAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `IDAtoms`}, optional
        if not `None`, then a list of `IDAtom` instance objects or an
        existing `IDAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return IDAtom

    def sort(self, key=attrgetter('mol', 'id'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def ids(self):
        """Return array of :attr:`IDAtom.id`\ s."""
        if len(set([atom.id for atom in self])) != len(self):
            self.assign_unique_ids()
        return np.asarray([atom.id for atom in self])

    @property
    def atom_ids(self):
        """Alias for :attr:`~IDAtoms.ids`."""
        return self.ids

    @property
    def serials(self):
        """Alias for :attr:`~IDAtoms.ids`."""
        return self.ids

    @property
    def mols(self):
        """Return array of `IDAtom.mol`\ s."""
        return np.asarray([atom.mol for atom in self])

    @property
    def mol_ids(self):
        """Alias for :attr:`~IDAtoms.mols`."""
        return self.mols

    @property
    def molecule_ids(self):
        """Alias for :attr:`~IDAtoms.mols`."""
        return self.mols

    @property
    def indices(self):
        """Return array of :attr:`IDAtom.index`\ s."""
        return np.asarray([self.index(atom) for atom in self])

    def assign_unique_ids(self, starting_id=1):
        """Assign unique :attr:`IDAtom.id` to each `IDAtom` in `IDAtoms`."""
        [setattr(atom, 'id', i) for i, atom in
         enumerate(self, start=starting_id)]

    def filter_ids(self, atom_ids, invert=False):
        """Filter `Atoms` by :attr:`IDAtoms.ids` in `atom_ids`.

        .. versionchanged:: 0.3.11

           Filters `Atoms` **in-place**. Use :meth:`~IDAtoms.filtered_ids`
           to get a new list of `Atoms`.

        Parameters
        ----------
        atom_ids : array_like
        invert : bool, optional

        """
        # mask = np.in1d(self.ids, atom_ids, invert=invert).nonzero()
        # self.data = np.asarray(self)[mask].tolist()
        self.data = [self.get_atom(id) for id in atom_ids]

    def filtered_ids(self, atom_ids, invert=False):
        """Return new `Atoms` object filtered by `atom_ids`.

        Returns a new `Atoms` object containing the `Atom` objects whose
        `Atom.id` is in `atom_ids` list.

        .. versionadded:: 0.3.11

        Parameters
        ----------
        atom_ids : array_like
        invert : bool, optional

        Returns
        -------
        filtered_atoms : `Atoms`
            An instance of `Atoms` (sub)class.

        """
        # mask = np.in1d(self.ids, atom_ids, invert=invert).nonzero()
        # return self.__class__(atoms=np.asarray(self)[mask].tolist(),
        #                       **self.kwargs)
        return self.__class__(atoms=[self.get_atom(id) for id in atom_ids],
                              **self.kwargs)

    def get_atom(self, id):
        """Get `IDAtom` with :attr:`Xatom.id` == `id`.

        Parameters
        ----------
        id : :class:`~python:int`

        Returns
        -------
        atom : `IDAtom` or `None`
            `IDAtom` instance if `IDAtoms` contains `IDAtom` with
            :attr:`IDAtom.id` == `id`, otherwise `None`

        """
        try:
            return self[np.where(self.ids == id)[0]]
        except TypeError:
            print('No atom with id = {}'.format(id))
            return None

    def get_atoms(self, ids=None, **kwargs):
        """Overrides parent class :meth:`Atoms.get_atoms`.

        Calls :meth:`~IDAtoms.filtered_ids` if `ids` is not None.

        Parameters
        ----------
        ids : {:class:`~python:list` of :class:`~python:int`\ s or \
            :class:`~python:None`}
        **kwargs : :class:`~python:dict`

        Returns
        -------
        :class:`~numpy:numpy.ndarray` if `asarray` is `True`
        :class:`~python:list` if `asarray` is `False` and `aslist` is `True`
        :class:`Atoms` object if `asarray` and `aslist` are `False`

        """
        if ids is None:
            return super().get_atoms(**kwargs)
        else:
            kwargs['aslist'] = kwargs.pop('aslist', False)
            return self.filtered_ids(ids).get_atoms(**kwargs)
