# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin classes for ring statistics (:mod:`sknano.core.atoms.mixins._ring_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins._ring_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import Counter, OrderedDict

import numpy as np
import pandas as pd

from sknano.core import timethis
# from sknano.core.analysis import _ring_finder
from sknano.core.analysis import find_rings
# from sknano.core.math import Vector

__all__ = ['RingAtomMixin', 'RingAtomsMixin']


class RingAtomMixin:
    """Mixin `Atom` class for ring statistics/network connectivity analysis.

    Attributes
    ----------
    Rn : :class:`~python:int`

    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.Rn = 0


class RingAtomsMixin:
    """Mixin `Atoms` class for ring statistics/network connectivity analysis.

    Attributes
    ----------
    ring_counter : :class:`~python:collections.Counter`
    rings : :class:`~python:collections.OrderedDict`

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

    @property
    def ring_stats(self):
        """:class:`~python:dict` of ring statistics."""
        try:
            return self._ring_stats
        except AttributeError:
            self.update_ring_stats()
            return self._ring_stats

    def update_ring_stats(self, angle_reference=None, bond_reference=None,
                          dihedral_reference=None, improper_reference=None):
        # stats = {}
        # stats['rings_per_atom'] = self.rings_per_atom
        data = []
        index = []
        rings = OrderedDict(sorted(self.rings.items()))
        for n, lor in rings.items():
            index.append(n)
            datadict = {}

            datadict['count'] = len(lor)
            datadict['centroids'] = centroids = []
            datadict['mean_angles'] = mean_angles = []
            datadict['mean_bonds'] = mean_bonds = []
            datadict['mean_angle_strains'] = mean_angle_strains = []
            datadict['mean_bond_strains'] = mean_bond_strains = []
            for ring in lor:
                centroids.append(ring.centroid)
                angles = ring.angles
                mean_angles.append(angles.mean)
                if angle_reference is not None:
                    angles.compute_strains(angle_reference)
                mean_angle_strains.append(np.mean(angles.strains))

                bonds = ring.bonds
                mean_bonds.append(bonds.mean)
                if bond_reference is not None:
                    bonds.compute_strains(bond_reference)
                mean_bond_strains.append(np.mean(bonds.strains))

                # dihedrals = ring.dihedrals
                # impropers = ring.impropers

            data.append(datadict)

        self._ring_stats = pd.DataFrame(data, index=index)

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
        self.reset_attrs(rings=True)
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

    def _update_ring_counter(self, ring_counts):
        """Private method for updating ring counter."""
        cntr = self.ring_counter
        for i, counts in enumerate(ring_counts):
            if counts > 0:
                cntr[i] = counts

    def _update_rings(self, nodes_list):
        for nodes in nodes_list:
            self.add_ring(nodes)

    def reset_attrs(self, rings=False, **kwargs):
        """Reset the :class:`RingAtomsMixin` class attributes, then call \
            parent class `reset_attrs` method."""
        if rings:
            self.reset_ring_atoms_attrs()
        super().reset_attrs(**kwargs)

    def reset_ring_atoms_attrs(self):
        """Reset the :class:`RingAtomsMixin` class attributes."""
        self.ring_counter = Counter()
        self.rings = OrderedDict()
        [setattr(atom, 'Rn', 0) for atom in self]

    def update_attrs(self, rings=False, **kwargs):
        """Update the :class:`RingAtomsMixin` class attributes."""
        super().update_attrs(**kwargs)
        if rings:
            self.analyze_network(**kwargs)
