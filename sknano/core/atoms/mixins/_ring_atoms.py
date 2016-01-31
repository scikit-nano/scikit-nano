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

from networkx.classes.graph import Graph

from sknano.core.math import Vector

__all__ = ['RingAtomMixin', 'RingAtomsMixin']


class RingAtomMixin:
    """Mixin `Atom` class for analysis of ring statistics/network \
        connectivity."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.visited = False
        self.rings = []
        self.linked_nodes = []
        self.Rc = 0
        self.Rn = 0

    # @property
    # def rings(self):
    #     return self._rings

    # @rings.setter
    # def rings(self, value):
    #     self._rings = value

    # @property
    # def visited(self):
    #     return self._visited

    # @visited.setter
    # def visited(self, value):
    #     if not isinstance(value, bool):
    #         raise ValueError('Expected boolean `True` or `False`')
    #     self._visited = value

    def visit(self, walk_vector=None):
        self.visited = True
        if walk_vector is None:
            walk_vector = Vector([0, 0, 0])
        # print('atom {} walk_vector: {}'.format(self.id, walk_vector))
        for nn in self.neighbors:
            if not nn.visited:
                nn.visit(walk_vector=walk_vector + (nn.r - self.r))


class RingAtomsMixin:
    """Mixin `Atoms` class for analysis of ring statistics/network \
        connectivity."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self.connectivity_matrix = None
        self.ring_counter = Counter()
        self.rings = None

    # @property
    # def linked_nodes(self):
    #     """2d :class:`~numpy:numpy.ndarray` of :attr:`Node.id`\ s."""
    #     return np.asarray([atom.linked_nodes for atom in self])

    # @property
    # def nlinked_nodes(self):
    #     return np.asarray([atom.nlinked_nodes for atom in self])

    @property
    def visited(self):
        return [atom.visited for atom in self]

    @property
    def graph(self):
        G = Graph()
        G.add_edges_from(self.bonds.indices)
        return G

    @property
    def rings_per_atom(self):
        return [atom.Rn for atom in self]

    def analyze_network(self, retcodes=None, **kwargs):
        """Analyze network properties.

        Parameters
        ----------
        retcodes : {None, :class:`~python:str`}, optional
        kwargs : :class:`~python:dict`

        Returns
        -------
        :class:`~python:tuple`
            if `retcodes` is None, the return value will be a tuple of
            :attr:`~RingAtomsMixin.ring_counter` and
            :attr:`~RingAtomsMixin.rings_per_atom`.

        """
        retvals = {}
        if retcodes is None:
            retcodes = ('ring_counter', 'rings_per_atom')
        rings = self.find_rings(cutoff=kwargs.get('cutoff', np.inf),
                                maxlength=kwargs.get('maxlength', 10))
        retvals['rings'] = rings
        retvals['ring_counter'] = self.ring_counter
        retvals['rings_per_atom'] = self.rings_per_atom
        return tuple(retvals[k] for k in retcodes)

    def find_rings(self, cutoff=np.inf, maxlength=10):
        """Find the shortest-path rings.

        Parameters
        ----------
        cutoff : :class:`~python:float`
        maxlength : :class:`~python:int`

        Returns
        -------
        :class:`~python:collections.Counter`

        """
        self.update_neighbor_lists()
        nn_idx = self.nn_idx
        seed = self.nn_seed_list
        r_nn = self.nn_vectors
        visited = np.zeros(self.NNN, dtype=bool)

        for i, atom in enumerate(self):
            for si in range(seed[i], seed[i+1]):
                ni = nn_idx[si]
                if i < ni:
                    visited[si] = True
                    for nsi in range(seed[ni], seed[ni+1]):
                        nni = nn_idx[nsi]
                        if nni == i:
                            visited[nsi] = True
                    self.traverse_network(root=i, node=ni, d=r_nn[si],
                                          visited=visited, maxlength=maxlength)

        return self.rings

    def traverse_network(self, root=None, node=None, d=None, visited=None,
                         maxlength=None, TOL=0.0001):
        """Traverse the network of atoms to determine its connectivity.

        Parameters
        ----------
        root, node : :class:`~python:int`
        d : array_like
        visited : array_like
        maxlength : :class:`~python:int`
        TOL : :class:`~python:float`

        """
        if maxlength is None:
            maxlength = -1

        seed = self.nn_seed_list
        amat = self.nn_adjacency_matrix
        r_nn = self.nn_vectors
        nn_idx = self.nn_idx
        cntr = self.ring_counter

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
                                        if amat[abs(nodes[m])][abs(nodes[n])] \
                                                != s:
                                            is_sp = False
                                if is_sp:
                                    for nidx in nodes:
                                        self[abs(nidx)].Rn += 1
                                    cntr[ring_size] += 1
                        elif dn == di - 1:
                            q.append((node, nodes + [-ni], d + r_nn[si]))
