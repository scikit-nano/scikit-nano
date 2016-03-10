# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
# distutils: language=c++

from libc.math cimport sqrt
from libc.stdlib cimport abs
from libcpp.deque cimport deque
from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

cdef struct queue_item:
    int32_t prev
    vector[int32_t] nodes
    float64_t r[3]


cdef void search_paths(int32_t root, int32_t node, int32_t node_si,
                       vector[int32_t] &visited, int32_t Natoms,
                       int32_t *nn_idx, int32_t *nn_seed, int32_t *nn_amat,
                       float64_t *nn_vecs, int32_t max_ring_size,
                       float64_t eps, vector[int32_t] &ring_counts,
                       vector[vector[int32_t]] &nodes_list):
    """Search for closed paths starting from `root` atom.

    Parameters
    ----------
    root, node : int32_t
    node_si : int32_t
    visited : int32_t array
    Natoms : int32_t
    nn_idx : int32_t array
    nn_seed : int32_t array
    nn_amat : int32_t array
    nn_vecs : float64_t array
    max_ring_size : int32_t
    eps : float64_t
    ring_counts : vector
    nodes_list : vector of vectors

    """
    cdef:
        int32_t i, j, ni, si, di, dni, n, m, s, abs_m, abs_n, ring_size
        int32_t prev, last
        float64_t r, rsq
        float64_t r0[3]
        bint found_ring, new_node = False
        vector[int32_t] nodes, tmpnodes
        queue_item qi
        deque[queue_item] q

    qi.prev = root
    qi.nodes.push_back(node)
    for i in range(3):
        qi.r[i] = nn_vecs[3 * node_si + i]
    q.push_back(qi)

    while q.size() > 0:
        qi = q.front()
        prev = qi.prev
        nodes = qi.nodes
        for i in range(3):
            r0[i] = qi.r[i]
        last = nodes.back()
        q.pop_front()
        if last > 0:
            i = last
            for si in range(nn_seed[i], nn_seed[i+1]):
                ni = nn_idx[si]
                if not visited[si] and ni != prev:
                    di = nn_amat[Natoms * root + i]
                    dni = nn_amat[Natoms * root + ni]
                    tmpnodes = nodes
                    if ((dni == di + 1) and
                        (nodes.size() < (max_ring_size - 1) / 2 or
                         max_ring_size < 0)):
                        tmpnodes.push_back(ni)
                        new_node = True
                    elif ((dni == di) or (dni == di - 1)):
                        tmpnodes.push_back(-ni)
                        new_node = True

                    if new_node:
                        qi.prev = last
                        qi.nodes = tmpnodes
                        for j in range(3):
                            qi.r[j] = r0[j] + nn_vecs[3 * si + j]
                        q.push_back(qi)
                        new_node = False
        else:
            i = -last
            for si in range(nn_seed[i], nn_seed[i + 1]):
                ni = nn_idx[si]
                if not visited[si] and ni != prev:
                    di = nn_amat[Natoms * root + i]
                    dni = nn_amat[Natoms * root + ni]
                    if ni == root:
                        rsq = 0.0
                        for j in range(3):
                            rsq += pow(r0[j] + nn_vecs[3 * si + j], 2)
                        r = sqrt(rsq)
                        # For a closed ring, the vector sum of the
                        # ring bond vectors will have zero length.
                        if r < eps:
                            found_ring = True
                            nodes.push_back(root)
                            ring_size = nodes.size()
                            for n in range(ring_size):
                                for m in range(n + 1, ring_size):
                                    s = m - n
                                    if s > ring_size / 2:
                                        s = ring_size - s

                                    abs_m = abs(nodes[m])
                                    abs_n = abs(nodes[n])

                                    if nn_amat[Natoms * abs_m + abs_n] != s:
                                        found_ring = False

                            if found_ring:
                                if ring_counts.size() < ring_size + 1:
                                    ring_counts.resize(ring_size + 1)
                                nodes_list.push_back(nodes)
                                ring_counts[ring_size] += 1

                    elif dni == di - 1:
                        tmpnodes = nodes
                        tmpnodes.push_back(-ni)

                        qi.prev = last
                        qi.nodes = tmpnodes
                        for j in range(3):
                            qi.r[j] = r0[j] + nn_vecs[3 * si + j]
                        q.push_back(qi)


cpdef tuple find_rings(int32_t Natoms, int32_t NNN, int32_t[::1] nn_idx,
                       int32_t[::1] nn_seed, int32_t[:,::1] nn_amat,
                       float64_t[:,::1] nn_vecs, int32_t max_ring_size,
                       float64_t eps):
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

    cdef:
        int32_t i, si, ni, nsi, nni
        vector[int32_t] visited
        vector[int32_t] ring_counts
        vector[vector[int32_t]] nodes_list

    visited.resize(NNN)

    for i in range(Natoms):
        for si in range(nn_seed[i], nn_seed[i+1]):
            ni = nn_idx[si]
            if i < ni:
                visited[si] = 1
                for nsi in range(nn_seed[ni], nn_seed[ni+1]):
                    nni = nn_idx[nsi]
                    if nni == i:
                        visited[nsi] = 1
                search_paths(i, ni, si, visited, Natoms, &nn_idx[0],
                             &nn_seed[0], &nn_amat[0,0], &nn_vecs[0,0],
                             max_ring_size, eps, ring_counts, nodes_list)

    return ring_counts, nodes_list
