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
from collections import Iterable, namedtuple
from functools import total_ordering
from operator import attrgetter

import numpy as np
# np.set_printoptions(edgeitems=20)
# np.set_printoptions(threshold=10000)
# import scipy as sp
from scipy import stats
# import pandas as pd

from sknano.core import BaseClass, UserList, TabulateMixin, rezero_array
from sknano.core.atoms import Atom, Atoms
# from sknano.core.math import vector as vec

__all__ = ['Topology', 'TopologyCollection', 'TopologyStats', 'check_operands',
           'AngularTopology', 'AngularTopologyCollection']


operand_error_msg = 'Expected an `iterable` object containing {}'
atoms_operand_error_msg = operand_error_msg.format('`Atom` objects')
ids_operand_error_msg = operand_error_msg.format('`ints`')
TopologyStats = namedtuple('TopologyStats', ('nobs', 'min', 'max', 'minmax',
                                             'mean', 'median', 'mode',
                                             'variance', 'std',
                                             'skewness', 'kurtosis'))


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
class Topology(BaseClass, TabulateMixin, metaclass=ABCMeta):
    """Base :class:`~sknano.core.atoms.Atom` topology class.

    Parameters
    ----------
    *atoms : {:class:`~python:list`, :class:`~sknano.core.atoms.Atoms`}
        :class:`~python:list` of :class:`~sknano.core.atoms.Atom`\ s
        or an :class:`~sknano.core.atoms.Atoms` object.
    size : :class:`~python:int`
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.
    id : :class:`~python:int`
    type : :class:`~python:int`
    check_operands : :class:`~python:bool`, optional

    Raises
    ------
    :class:`~python:TypeError`
        if `atoms` is not a list of :class:`~sknano.core.atoms.Atom` objects
        or an :class:`~sknano.core.atoms.Atoms` object.
    :class:`~python:ValueError`
        if len(atoms) != `size`.
    """
    def __init__(self, *atoms, size, id=0, type=0, parent=None,
                 check_operands=True):
        if check_operands:
            if not isinstance(atoms, Iterable):
                raise TypeError(atoms_operand_error_msg)
            if len(atoms) == 1:
                atoms = atoms[0]
            if not isinstance(atoms, (Iterable, Atoms)) or not \
                    all([isinstance(atom, Atom) for atom in atoms]):
                raise TypeError(atoms_operand_error_msg)
            if len(atoms) != size:
                raise ValueError('Expected {} atoms'.format(size))

        super().__init__()
        from .. import StructureAtoms

        self.atoms = StructureAtoms()
        self.atoms.extend(list(atoms))

        self.check_operands = check_operands
        self.id = id
        self.type = type
        self.parent = parent
        self.fmtstr = "id={id!r}, type={type!r}, parent={parent!r}"

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
        return ['atoms', 'measure', 'id', 'type', 'parent']

    # def __iter__(self):
    #     return iter(self.atoms)

    @property
    def atoms(self):
        """:class:`~sknano.core.atoms.Atoms` in `TopologyCollection`."""
        return self._atoms

    @atoms.setter
    def atoms(self, value):
        # TODO: This should perform some sort of type check as well as update
        # other attributes.
        self._atoms = value

    @property
    def atom_ids(self):
        """:attr:`Topology.atoms` :attr:`~sknano.core.atoms.IDAtoms.ids`."""
        return tuple(self.atoms.ids)

    @property
    def centroid(self):
        """:attr:`~sknano.core.atoms.XYZAtoms.centroid` of \
            :attr:`Topology.atoms`."""
        return self.atoms.centroid

    @property
    def measure(self):
        """Measure of topology."""
        try:
            return self._measure
        except AttributeError:
            self._update_measure()
            return self._measure

    @property
    def strain(self):
        """Strain in measure."""
        try:
            return self._strain
        except AttributeError:
            return 0.0

    @abstractmethod
    def compute_measure(self):
        """Compute topological measure from :attr:`Topology.atoms`."""
        raise NotImplementedError

    def compute_strain(self, m0):
        """Compute topological strain in :attr:`Topology.measure`.

        Parameters
        ----------
        m0 : :class:`~python:float`

        Returns
        -------
        :class:`~python:float`

        """
        m = self.measure
        self._strain = (m0 - m) / m
        return self._strain

    def rotate(self, **kwargs):
        """Rotate the `Topology` by rotating the \
            :attr:`~Topology.atoms`."""
        [atom.rotate(fix_anchor_point=True, **kwargs) for atom in self.atoms]

    def _update_measure(self):
        self._measure = self.compute_measure()

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(id=self.id, type=self.type, parent=self.parent)


class TopologyCollection(TabulateMixin, UserList):
    """Base :class:`~sknano.core.atoms.Atoms` topology class.

    Parameters
    ----------
    topolist : {None, sequence, `Topology`}, optional
        if not `None`, then a list of `Topology` objects
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.

    """
    def __init__(self, topolist=None, parent=None):
        super().__init__(initlist=topolist)
        self.parent = parent
        self.fmtstr = super().fmtstr + ", parent={parent!r}"

    @property
    def __item_class__(self):
        return Topology

    def sort(self, key=attrgetter('measure'), reverse=False):
        super().sort(key=key, reverse=reverse)

    @property
    def Ntypes(self):
        """Number of unique :attr:`Topology.type`\ s."""
        return len(set(self.types))

    @property
    def parent(self):
        """Parent :class:`~sknano.core.atoms.Molecule`, if any."""
        return self._parent

    @parent.setter
    def parent(self, value):
        self._parent = self.kwargs['parent'] = value

    # @property
    # def atoms(self):
    #     """`Atoms` :class:`python:set` in `TopologyCollection`."""
    #     atoms = []
    #     [atoms.extend(topology.atoms) for topology in self]
    #     atoms = \
    #         list(dedupe(list(flatten([topology.atoms for topology in self])),
    #                     key=attrgetter('id')))
    #     return StructureAtoms(atoms)

    @property
    def measures(self):
        """:class:`~numpy:numpy.ndarray` of \
            :attr:`~Topology.measure`\ s."""
        try:
            return self._measures
        except AttributeError:
            self._update_measures()
            return self._measures

    @property
    def strains(self):
        """:class:`~numpy:numpy.ndarray` of \
            :attr:`~Topology.strain`\ s."""
        try:
            return self._strains
        except AttributeError:
            return np.zeros(len(self), dtype=float)

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
        """An alias for :attr:`TopologyCollection.mean_measure`."""
        return self.mean_measure

    @property
    def atom_ids(self):
        """:class:`~python:list` of :attr:`~Topology.atom_ids`."""
        # return np.asarray([topology.atom_ids for topology in self])
        return [topology.atom_ids for topology in self]

    @property
    def ids(self):
        """:class:`~python:list` of :attr:`~Topology.id`\ s."""
        # return np.asarray([topology.id for topology in self])
        return [topology.id for topology in self]

    @property
    def types(self):
        """:class:`~python:list` of :attr:`~Topology.type`\ s."""
        # return np.asarray([topology.type for topology in self])
        return [topology.type for topology in self]

    @property
    def unique(self):
        """Return new `TopologyCollection` object containing the set of \
            unique `Topology`\ s."""
        seen = set()
        unique = []
        for topology in self:
            atom_ids = tuple(topology.atom_ids)
            ratom_ids = tuple(reversed(atom_ids))
            if atom_ids not in seen and ratom_ids not in seen:
                unique.append(topology)
                seen.add(atom_ids)
                seen.add(ratom_ids)
        return self.__class__(topolist=unique, **self.kwargs)

    @property
    def statistics(self):
        """:class:`TopologyStats` of :attr:`TopologyCollection.measures` \
            statistics."""
        # measures = self.unique.measures
        measures = self.measures
        topostats = stats.describe(measures)._asdict()
        topostats['min'] = np.min(measures)
        topostats['max'] = np.max(measures)
        topostats['median'] = np.median(measures)
        topostats['mode'] = stats.mode(measures).mode[0]
        topostats['std'] = np.std(measures)
        return TopologyStats(**topostats)

    def compute_strains(self, reference):
        """Return :class:`~numpy:numpy.ndarray` of topology measure strains.

        Parameters
        ----------
        reference : :class:`~python:float` or array_like
            Reference/starting measure.

        Returns
        -------
        :class:`~numpy:numpy.ndarray`

        """
        if not (np.isscalar(reference) or
                isinstance(reference, (list, np.ndarray))):
            raise TypeError('Expected a `float` or `array_like` object.')
        if np.isscalar(reference):
            self._strains = \
                np.asarray([topology.compute_strain(reference) for
                            topology in self])
        else:
            if len(reference) != len(self):
                raise ValueError('`reference` must be same length as '
                                 'topology list')
            self._strains = \
                np.asarray([topology.compute_strain(ref) for topology, ref in
                            zip(self, reference)])
        return self._strains

    def _update_measures(self):
        [topology._update_measure() for topology in self]
        self._measures = \
            rezero_array(np.asarray([topology.measure for topology in self]))

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        super_dict = super().todict()
        super_dict.update(dict(parent=self.parent))
        return super_dict


class AngularTopology(Topology):
    """`Topology` sub-class for topology with angular measure.

    Parameters
    ----------
    *atoms : {:class:`~python:list`, :class:`~sknano.core.atoms.Atoms`}
        :class:`~python:list` of :class:`~sknano.core.atoms.Atom`\ s
        or an :class:`~sknano.core.atoms.Atoms` object.
    check_operands : :class:`~python:bool`, optional
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.
    id : :class:`~python:int`
    degrees : :class:`~python:bool`, optional

    Raises
    ------
    :class:`~python:TypeError`
        if `atoms` is not a list of :class:`~sknano.core.atoms.Atom` objects
        or an :class:`~sknano.core.atoms.Atoms` object.
    """
    def __init__(self, *args, degrees=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.degrees = degrees
        self.fmtstr = super().fmtstr + ", degrees={degrees!r}"

    @property
    def angle(self):
        """An alias for :attr:`Topology.measure`."""
        return self.measure

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        super_dict = super().todict()
        super_dict.update(dict(degrees=self.degrees))
        return super_dict


class AngularTopologyCollection(TopologyCollection):
    """`TopologyCollection` sub-class for collection of angular topologies.

    Parameters
    ----------
    topolist : {None, sequence, `AngularTopology`}, optional
        if not `None`, then a list of `AngularTopology` objects
    parent : Parent :class:`~sknano.core.atoms.Molecule`, if any.
    degrees : :class:`~python:bool`, optional

    """
    def __init__(self, *args, degrees=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.degrees = degrees
        self.fmtstr = super().fmtstr + ", degrees={degrees!r}"

    @property
    def __item_class__(self):
        return AngularTopology

    @property
    def degrees(self):
        """:class:`~python:bool` setting for returning angles in degrees."""
        return self._degrees

    @degrees.setter
    def degrees(self, value):
        if not isinstance(value, bool):
            raise ValueError('Expected a boolean value.')
        self._degrees = self.kwargs['degrees'] = value
        [setattr(topoobj, 'degrees', value) for topoobj in self]
        super()._update_measures()

    @property
    def angles(self):
        """:class:`~numpy:numpy.ndarray` of \
            :attr:`~AngularTopology.angle`\ s."""
        return self.measures

    @property
    def mean_angle(self):
        """Mean angle."""
        return self.mean_measure

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        super_dict = super().todict()
        super_dict.update(dict(degrees=self.degrees))
        return super_dict
