# -*- coding: utf-8 -*-
"""
===============================================================================
Ring analysis (:mod:`sknano.core.analysis.ring_analysis`)
===============================================================================

.. currentmodule:: sknano.core.analysis.ring_analysis

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import deque

import numpy as np

# from sknano.core import timethis
# from sknano.core.math import Vector

try:
    from . import _ring_finder
except ImportError:
    _ring_finder = None

__all__ = ['find_rings', 'py_find_rings']


# @timethis
def find_rings(Natoms, NNN, nn_idx, nn_seed, nn_amat, nn_vecs,
               max_ring_size=None, eps=0.0001, pyversion=False):
    """Find rings.

    Parameters
    ----------
    Natoms : :class:`~python:int`
    NNN : :class:`~python:int`
    nn_idx : :class:`~numpy:numpy.ndarray`
    nn_seed : :class:`~numpy:numpy.ndarray`
    nn_amat : :class:`~numpy:numpy.ndarray`
    nn_vecs : :class:`~numpy:numpy.ndarray`
    max_ring_size : :class:`~python:int`, optional
    eps : :class:`~python:float`, optional

    """
    if max_ring_size is None:
        max_ring_size = -1

    if _ring_finder is None or pyversion:
        return py_find_rings(Natoms, NNN, nn_idx, nn_seed, nn_amat, nn_vecs,
                             max_ring_size, eps)

    nn_idx = np.asarray(nn_idx, dtype=np.intc)
    nn_seed = np.asarray(nn_seed, dtype=np.intc)
    nn_amat = np.asarray(nn_amat, dtype=np.intc)
    try:
        nn_vecs = np.asarray([v.__array__() for v in nn_vecs], dtype=np.double)
    except AttributeError:
        nn_vecs = np.asarray(nn_vecs, dtype=np.double)

    return _ring_finder.find_rings(Natoms, NNN, nn_idx, nn_seed,
                                   nn_amat, nn_vecs, max_ring_size, eps)


# @timethis
def py_find_rings(Natoms, NNN, nn_idx, nn_seed, nn_amat, nn_vecs,
                  max_ring_size=None, eps=0.0001):
    """Find rings.

    Parameters
    ----------
    Natoms : :class:`~python:int`
    NNN : :class:`~python:int`
    nn_idx : :class:`~numpy:numpy.ndarray`
    nn_seed : :class:`~numpy:numpy.ndarray`
    nn_amat : :class:`~numpy:numpy.ndarray`
    nn_vecs : :class:`~numpy:numpy.ndarray`
    max_ring_size : :class:`~python:int`, optional
    eps : :class:`~python:float`, optional

    """
    if max_ring_size is None:
        max_ring_size = -1

    visited = np.zeros(NNN, dtype=bool)
    # ring_counter = Counter()
    ring_counts = []
    nodes_list = []

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
                for si in range(nn_seed[i], nn_seed[i+1]):
                    ni = nn_idx[si]
                    if not visited[si] and ni != prev:
                        di = nn_amat[root][i]
                        dni = nn_amat[root][ni]
                        if (dni == di + 1):
                            if (len(nodes) < (max_ring_size - 1) // 2 or
                                    max_ring_size < 0):
                                q.append((node, nodes + [ni],
                                          d + nn_vecs[si]))
                            else:
                                continue
                        elif ((dni == di) or (dni == di - 1)):
                            q.append((node, nodes + [-ni], d + nn_vecs[si]))
                        else:
                            raise ValueError('bad adjacency matrix')
            else:
                i = -node
                for si in range(nn_seed[i], nn_seed[i + 1]):
                    ni = nn_idx[si]
                    if not visited[si] and ni != prev:
                        di = nn_amat[root][i]
                        dni = nn_amat[root][ni]
                        if ni == root:
                            dr = nn_vecs[si] + d
                            # For a closed ring, the vector sum of the
                            # ring bond vectors will have zero length.
                            if dr.norm < eps:
                                found_ring = True
                                nodes.append(root)
                                ring_size = len(nodes)
                                for n in range(ring_size):
                                    for m in range(n + 1, ring_size):
                                        s = m - n
                                        if s > ring_size // 2:
                                            s = ring_size - s
                                        if nn_amat[abs(nodes[m])][
                                                abs(nodes[n])] != s:
                                            found_ring = False
                                if found_ring:
                                    nodes_list.append(nodes)
                                    while len(ring_counts) < ring_size + 1:
                                        ring_counts.append(0)
                                    ring_counts[ring_size] += 1
                        elif dni == di - 1:
                            q.append((node, nodes + [-ni], d + nn_vecs[si]))

    for i in range(Natoms):
        for si in range(nn_seed[i], nn_seed[i+1]):
            ni = nn_idx[si]
            if i < ni:
                visited[si] = True
                for nsi in range(nn_seed[ni], nn_seed[ni+1]):
                    nni = nn_idx[nsi]
                    if nni == i:
                        visited[nsi] = True
                search_paths(root=i, node=ni, d=nn_vecs[si])

    return ring_counts, nodes_list
