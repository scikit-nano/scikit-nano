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

from collections import Counter, deque

import numpy as np

# from sknano.core.math import Vector

__all__ = ['RingAtomMixin', 'RingAtomsMixin']


class RingAtomMixin:
    """Mixin `Atom` class for analysis of ring statistics/network \
        connectivity."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.rings = []
        self.linked_nodes = []
        self.Rn = 0


class RingAtomsMixin:
    """Mixin `Atoms` class for analysis of ring statistics/network \
        connectivity."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ring_counter = Counter()
        self.rings = None

    @property
    def rings_per_atom(self):
        return [atom.Rn for atom in self]

    @property
    def Rn_counter(self):
        """Return rings per atom counter.

        Returns
        -------
        :class:`~python:collections.Counter`

        """
        return Counter([atom.Rn for atom in self])

    def analyze_network(self, cutoff=np.inf, maxlength=None, TOL=0.0001,
                        retcodes=None):
        """Find the shortest-path rings.

        Parameters
        ----------
        cutoff : :class:`~python:float`, optional
        maxlength : :class:`~python:int`, optional
        TOL : :class:`~python:float`, optional
        retcodes : :class:`~python:tuple`, optional

        """
        if maxlength is None:
            maxlength = -1

        self.update_neighbor_lists()
        nn_idx = self.nn_idx
        seed = self.nn_seed_list
        r_nn = self.nn_vectors
        visited = np.zeros(self.NNN, dtype=bool)
        amat = self.nn_adjacency_matrix
        cntr = self.ring_counter

        def search_paths(root=None, node=None, d=None):
            """Search for closed paths starting from `root` atom.

            Parameters
            ----------
            root, node : :class:`~python:int`
            d : array_like

            """
            q = deque([(root, [node], d)])
            while q:
                prev, nodes, d = q.popleft()
                node = nodes[-1]
                if node > 0:
                    i = node
                    for si in range(seed[i], seed[i+1]):
                        ni = nn_idx[si]
                        if not visited[si] and ni != prev:
                            di = amat[root][i]
                            dn = amat[root][ni]
                            if ((dn == di + 1) and
                                (len(nodes) < (maxlength - 1) / 2
                                 or maxlength < 0)):
                                q.append((node, nodes + [ni], d + r_nn[si]))
                            elif ((dn == di) or (dn == di - 1)):
                                q.append((node, nodes + [-ni], d + r_nn[si]))
                else:
                    i = -node
                    for si in range(seed[i], seed[i + 1]):
                        ni = nn_idx[si]
                        if not visited[si] and ni != prev:
                            di = amat[root][i]
                            dn = amat[root][ni]
                            if ni == root:
                                dr = r_nn[si] + d
                                if dr.norm < TOL:
                                    is_sp = True
                                    nodes.append(root)
                                    ring_size = len(nodes)
                                    for n in range(ring_size):
                                        for m in range(n + 1, ring_size):
                                            s = m - n
                                            if s > ring_size / 2:
                                                s = ring_size - s
                                            if amat[abs(nodes[m])][
                                                    abs(nodes[n])] != s:
                                                is_sp = False
                                    if is_sp:
                                        for nidx in nodes:
                                            self[abs(nidx)].Rn += 1
                                        cntr[ring_size] += 1
                            elif dn == di - 1:
                                q.append((node, nodes + [-ni], d + r_nn[si]))

        for i, atom in enumerate(self):
            for si in range(seed[i], seed[i+1]):
                ni = nn_idx[si]
                if i < ni:
                    visited[si] = True
                    for nsi in range(seed[ni], seed[ni+1]):
                        nni = nn_idx[nsi]
                        if nni == i:
                            visited[nsi] = True
                    search_paths(root=i, node=ni, d=r_nn[si])

        if retcodes is not None and isinstance(retcodes, (tuple, list)):
            return tuple([getattr(self, attr, None) for attr in retcodes])
