# -*- coding: utf-8 -*-
"""
===============================================================================
Ring analysis (:mod:`sknano.core.analysis._ring_analysis`)
===============================================================================

.. currentmodule:: sknano.core.analysis._ring_analysis

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import Counter, deque

import numpy as np

# from sknano.core import timethis
# from sknano.core.math import Vector

try:
    from . import _ring_finder
except ImportError:
    _ring_finder = None

__all__ = ['find_rings', 'py_find_rings']


def find_rings(Natoms, NNN, nn_idx, seed, amat, r_nn, maxlength=None,
               eps=0.0001, pyversion=False):
    """Analyze the network connectivity.

    Parameters
    ----------
    cutoff : :class:`~python:float`, optional
    maxlength : :class:`~python:int`, optional
    eps : :class:`~python:float`, optional

    """
    if maxlength is None:
        maxlength = -1

    if _ring_finder is None or pyversion:
        return py_find_rings(Natoms, NNN, nn_idx, seed, amat, r_nn,
                             maxlength, eps)

    nn_idx = np.asarray(nn_idx, dtype=np.intc)
    seed = np.asarray(seed, dtype=np.intc)
    amat = np.asarray(amat, dtype=np.intc)
    try:
        r_nn = np.asarray([v.__array__() for v in r_nn], dtype=np.double)
    except AttributeError:
        r_nn = np.asarray(r_nn, dtype=np.double)

    return _ring_finder.find_rings(Natoms, NNN, nn_idx, seed,
                                   amat, r_nn, maxlength, eps)


def py_find_rings(Natoms, NNN, nn_idx, seed, amat, r_nn, maxlength=None,
                  eps=0.0001):
    """Analyze the network connectivity.

    Parameters
    ----------
    cutoff : :class:`~python:float`, optional
    maxlength : :class:`~python:int`, optional
    eps : :class:`~python:float`, optional

    """
    if maxlength is None:
        maxlength = -1

    visited = np.zeros(NNN, dtype=bool)
    ring_counter = Counter()
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
                for si in range(seed[i], seed[i+1]):
                    ni = nn_idx[si]
                    if not visited[si] and ni != prev:
                        di = amat[root][i]
                        dni = amat[root][ni]
                        if (dni == di + 1):
                            if (len(nodes) < (maxlength - 1) // 2 or
                                    maxlength < 0):
                                q.append((node, nodes + [ni],
                                          d + r_nn[si]))
                            else:
                                continue
                        elif ((dni == di) or (dni == di - 1)):
                            q.append((node, nodes + [-ni], d + r_nn[si]))
                        else:
                            raise ValueError('bad adjacency matrix')
            else:
                i = -node
                for si in range(seed[i], seed[i + 1]):
                    ni = nn_idx[si]
                    if not visited[si] and ni != prev:
                        di = amat[root][i]
                        dni = amat[root][ni]
                        if ni == root:
                            dr = r_nn[si] + d
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
                                        if amat[abs(nodes[m])][
                                                abs(nodes[n])] != s:
                                            found_ring = False
                                if found_ring:
                                    nodes_list.append(nodes)
                                    ring_counter[ring_size] += 1
                        elif dni == di - 1:
                            q.append((node, nodes + [-ni], d + r_nn[si]))

    for i in range(Natoms):
        for si in range(seed[i], seed[i+1]):
            ni = nn_idx[si]
            if i < ni:
                visited[si] = True
                for nsi in range(seed[ni], seed[ni+1]):
                    nni = nn_idx[nsi]
                    if nni == i:
                        visited[nsi] = True
                search_paths(root=i, node=ni, d=r_nn[si])

    return ring_counter, nodes_list
