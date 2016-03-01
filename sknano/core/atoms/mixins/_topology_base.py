# -*- coding: utf-8 -*-
"""
===============================================================================
Base topology classes (:mod:`sknano.core.atoms.mixins._topology_base`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins._topology_base

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numbers
from abc import ABCMeta, abstractmethod
from collections import Iterable
from functools import total_ordering
from operator import attrgetter

import numpy as np
# np.set_printoptions(edgeitems=20)
np.set_printoptions(threshold=10000)

from sknano.core import BaseClass, UserList, rezero_array
from sknano.core.atoms import Atom, Atoms
# from sknano.core.math import vector as vec

__all__ = ['AtomTopology', 'AtomsTopology', 'check_operands']


operand_error_msg = 'Expected an `iterable` object containing {}'
atoms_operand_error_msg = operand_error_msg.format('`Atom` objects')
ids_operand_error_msg = operand_error_msg.format('`ints`')


def check_operands(*atoms, size=None):
    """Check atom operands.

    Parameters
    ----------
    *atoms : {:class:`~python:list`, :class:`~sknano.core.atoms.Atoms`}
        :class:`~python:list` of :class:`~sknano.core.atoms.Atom`\ s
        or an :class:`~sknano.core.atoms.Atoms` object.
    size : :class:`~python:int`

    Returns
    -------
    :class:`~python:tuple`

    Raises
    ------
    :class:`~python:TypeError`
        if `atoms` is not a list of :class:`~sknano.core.atoms.Atom` objects
        or an :class:`~sknano.core.atoms.Atoms` object.
    :class:`~python:ValueError`
        if len(atoms) != `size`.

    """
    if size is None:
        raise ValueError('Expected `int` for size')
    if not isinstance(atoms, Iterable):
        raise TypeError(atoms_operand_error_msg)
    if len(atoms) == 1:
        if isinstance(atoms[0], Iterable):
            return check_operands(*atoms[0], size=size)
        else:
            atoms = atoms[0]
    if not isinstance(atoms, (Iterable, Atoms)) or not \
            all([isinstance(atom, Atom) for atom in atoms]):
        raise TypeError(atoms_operand_error_msg)
    if len(atoms) != size:
        raise ValueError('Expected {} atoms'.format(size))

    return atoms


@total_ordering
class AtomTopology(BaseClass, metaclass=ABCMeta):
    """Base :class:`~sknano.core.atoms.Atom` topology class.

    Parameters
    ----------
    *atoms : {:class:`~python:list`, :class:`~sknano.core.atoms.Atoms`}
        :class:`~python:list` of :class:`~sknano.core.atoms.Atom`\ s
        or an :class:`~sknano.core.atoms.Atoms` object.
    check_operands : :class:`~python:bool`, optional
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.
    id : :class:`~python:int`

    Raises
    ------
    :class:`~python:TypeError`
        if `atoms` is not a list of :class:`~sknano.core.atoms.Atom` objects
        or an :class:`~sknano.core.atoms.Atoms` object.
    """
    def __init__(self, *atoms, check_operands=True, parent=None, id=0):
        if check_operands:
            if not isinstance(atoms, Iterable):
                raise TypeError(atoms_operand_error_msg)
            if len(atoms) == 1:
                atoms = atoms[0]
            if not isinstance(atoms, (Iterable, Atoms)) or not \
                    all([isinstance(atom, Atom) for atom in atoms]):
                raise TypeError(atoms_operand_error_msg)

        super().__init__()
        from .. import StructureAtoms

        self.atoms = StructureAtoms()
        self.atoms.extend(list(atoms))

        self.check_operands = check_operands
        self.parent = parent
        self.id = id
        self.fmtstr = "parent={parent!r}, id={id!r}"

    def _is_valid_operand(self, other):
        return isinstance(other, (numbers.Number, self.__class__))

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if np.isscalar(other):
            return np.allclose(self.measure, other)
        return self is other or self.atoms == other.atoms

    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if np.isscalar(other):
            return self.measure < other
        return self.atoms < other.atoms

    def __dir__(self):
        return ['atoms', 'measure', 'parent', 'id']

    # def __iter__(self):
    #     return iter(self.atoms)

    @property
    def atoms(self):
        """:class:`~sknano.core.atoms.Atoms` in `AtomsTopology`."""
        return self._atoms

    @atoms.setter
    def atoms(self, value):
        # TODO: This should perform some sort of type check as well as update
        # other attributes.
        self._atoms = value

    @property
    def ids(self):
        """:attr:`AtomTopology.atoms` \
            :attr:`~sknano.core.atoms.IDAtoms.ids`."""
        return tuple(self.atoms.ids)

    @property
    def centroid(self):
        """:attr:`~sknano.core.atoms.XYZAtoms.centroid` of \
            :attr:`AtomTopology.atoms`."""
        return self.atoms.centroid

    @property
    def measure(self):
        """Measure of topology."""
        return self.compute_measure()

    @abstractmethod
    def compute_measure(self):
        """Compute topological measure from :attr:`AtomTopology.atoms`."""
        raise NotImplementedError

    def compute_strain(self, m0):
        """Compute topological strain in :attr:`AtomTopology.measure`.

        Parameters
        ----------
        m0 : :class:`~python:float`

        Returns
        -------
        :class:`~python:float`

        """
        m = self.measure
        return (m0 - m) / m

    def rotate(self, **kwargs):
        """Rotate the `AtomTopology` by rotating the \
            :attr:`~AtomTopology.atoms`."""
        [atom.rotate(fix_anchor_point=True, **kwargs) for atom in self.atoms]

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(parent=self.parent, id=self.id)


class AtomsTopology(UserList):
    """Base :class:`~sknano.core.atoms.Atoms` topology class.

    Parameters
    ----------
    topolist : {None, sequence, `AtomTopology`}, optional
        if not `None`, then a list of `AtomTopology` objects
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.

    """
    def __init__(self, topolist=None, parent=None):
        super().__init__(initlist=topolist)
        self.parent = parent
        self.fmtstr = super().fmtstr + ", parent={parent!r}"

    @property
    def __item_class__(self):
        return AtomTopology

    def sort(self, key=attrgetter('measure'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def parent(self):
        """Parent :class:`~sknano.core.atoms.Molecule`, if any."""
        return self._parent

    @parent.setter
    def parent(self, value):
        self._parent = self.kwargs['parent'] = value

    # @property
    # def atoms(self):
    #     """`Atoms` :class:`python:set` in `AtomsTopology`."""
    #     atoms = []
    #     [atoms.extend(topology.atoms) for topology in self]
    #     atoms = \
    #         list(dedupe(list(flatten([topology.atoms for topology in self])),
    #                     key=attrgetter('id')))
    #     return StructureAtoms(atoms)

    @property
    def measures(self):
        """:class:`~numpy:numpy.ndarray` of \
            :attr:`~AtomTopology.measure`\ s."""
        return rezero_array(
            np.asarray([topology.measure for topology in self]))

    @property
    def mean_measure(self):
        """Mean measure."""
        if np.all(self.measures == np.inf):
            return np.inf
        if np.any(self.measures == np.inf):
            return np.ma.mean(np.ma.array(self.measures,
                                          mask=self.measures == np.inf))
        return np.mean(self.measures)

    @property
    def mean(self):
        """An alias for :attr:`AtomsTopology.mean_measure`."""
        return self.mean_measure

    @property
    def ids(self):
        """Return array of :attr:`~AtomTopology.ids`."""
        # return np.asarray([topology.ids for topology in self])
        return [topology.ids for topology in self]

    @property
    def unique(self):
        """Return new `AtomsTopology` object containing the set of unique \
            `AtomTopology`\ s."""
        seen = set()
        unique = []
        for topology in self:
            if topology.ids not in seen and \
                    tuple(reversed(topology.ids)) not in seen:
                unique.append(topology)
                seen.add(topology.ids)
                seen.add(tuple(reversed(topology.ids)))
        return self.__class__(topolist=unique)

    @property
    def statistics(self):
        """:class:`~python:dict` of :class:`AtomsTopology.measures` \
            statistics."""
        stats = {}
        measures = self.unique.measures
        stats['min'] = np.min(measures)
        stats['max'] = np.max(measures)
        stats['mean'] = np.mean(measures)
        return stats

    def compute_strain(self, m0):
        """Return :class:`~numpy:numpy.ndarray` of topology measure strains.

        Parameters
        ----------
        m0 : :class:`~python:float`
            Reference/starting measure.

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        """
        return np.asarray([topology.compute_strain(m0) for topology in self])

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        super_dict = super().todict()
        super_dict.update(dict(parent=self.parent))
        return super_dict
