# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes for ring statistics (:mod:`sknano.core.atoms.mixins._ring_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins._ring_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import Counter, OrderedDict

import numpy as np

from sknano.core import timethis
# from sknano.core.math import Vector

# from sknano.core.analysis import _ring_finder
from sknano.core.analysis import find_rings


__all__ = ['RingAtomMixin', 'RingAtomsMixin']


class RingAtomMixin:
    """Mixin `Atom` class for ring statistics/network connectivity analysis.

    Attributes
    ----------
    ring_counter : :class:`~python:collections.Counter`
    rings : :class:`~python:list`

    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.rings = []
        self.linked_nodes = []
        self.Rn = 0


class RingAtomsMixin:
    """Mixin `Atoms` class for ring statistics/network connectivity analysis.

    Attributes
    ----------
    ring_counter : :class:`~python:collections.Counter`
    rings : :class:`~python:list`

    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ring_counter = Counter()
        self.rings = OrderedDict()

    @property
    def rings_per_atom(self):
        """Rings per atom list."""
        return [atom.Rn for atom in self]

    @property
    def Rn_counter(self):
        """Return rings per atom counter.

        Returns
        -------
        :class:`~python:collections.Counter`

        """
        return Counter([atom.Rn for atom in self])

    def add_ring(self, nodes):
        """Append ring atoms to :attr:`~RingAtomsMixin.rings`."""
        for nidx in nodes:
            self[abs(nidx)].Rn += 1
        n = len(nodes)
        rings = self.rings
        if n not in rings:
            rings[n] = []

        rings[n].append(self.__class__([self[abs(nidx)] for nidx in nodes],
                                       update_item_class=False, **self.kwargs))

    def _update_ring_counter(self, ring_counts):
        cntr = self.ring_counter
        for i, counts in enumerate(ring_counts):
            if counts > 0:
                cntr[i] = counts

    def _update_rings(self, nodes_list):
        for nodes in nodes_list:
            self.add_ring(nodes)

    def reset_attrs(self):
        """Reset the :class:`~RingAtomsMixin` class attributes then call \
            super class method."""
        self.reset_ring_atoms_attrs()
        super().reset_attrs()

    def reset_ring_atoms_attrs(self):
        """Reset the :class:`~RingAtomsMixin` class attributes."""
        self.ring_counter = Counter()
        self.rings = OrderedDict()
        [setattr(atom, 'Rn', 0) for atom in self]

    @timethis
    def analyze_network(self, cutoff=np.inf, max_ring_size=None, eps=0.0001,
                        retcodes=None, pyversion=False):
        """Analyze the network connectivity.

        Parameters
        ----------
        cutoff : :class:`~python:float`, optional
        max_ring_size : :class:`~python:int`, optional
        eps : :class:`~python:float`, optional
        retcodes : :class:`~python:tuple`, optional

        """
        if max_ring_size is None:
            max_ring_size = -1
        self.reset_ring_atoms_attrs()
        if self.verbose:
            print('Analyzing the network connectivity...')

        if not self.neighbors_analyzed:
            self.update_neighbors(cutoffs=[cutoff])

        self.update_neighbor_lists()
        Natoms = self.Natoms
        NNN = self.NNN
        nn_idx = self.nn_idx
        nn_seed = self.nn_seed
        nn_amat = self.nn_adjacency_matrix
        nn_vecs = self.nn_vectors

        ring_counts, nodes_list = \
            find_rings(Natoms, NNN, nn_idx, nn_seed, nn_amat, nn_vecs,
                       max_ring_size, eps, pyversion)

        self._update_ring_counter(ring_counts)
        self._update_rings(nodes_list)

        if retcodes is not None and isinstance(retcodes, (tuple, list)):
            return tuple([getattr(self, attr, None) for attr in retcodes])
