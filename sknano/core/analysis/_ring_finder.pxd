from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

ctypedef np.float32_t float32_t
ctypedef np.float64_t float64_t
ctypedef np.int32_t int32_t
ctypedef np.int64_t int64_t
ctypedef np.uint32_t uint32_t
ctypedef np.uint64_t uint64_t

cdef struct queue_item


cdef void search_paths(int32_t root, int32_t node, int32_t node_si,
                       vector[int32_t] &visited, int32_t Natoms,
                       int32_t *nn_idx, int32_t *nn_seed, int32_t *nn_amat,
                       float64_t *nn_vecs, int32_t max_ring_size, float64_t eps,
                       vector[int32_t] &ring_counts,
                       vector[vector[int32_t]] &nodes_list)


cpdef tuple find_rings(int32_t Natoms, int32_t NNN, int32_t[::1] nn_idx,
                       int32_t[::1] nn_seed, int32_t[:,::1] nn_amat,
                       float64_t[:,::1] nn_vecs, int32_t max_ring_size,
                       float64_t eps)
