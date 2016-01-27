# -*- coding: utf-8 -*-
"""
==============================================================================
Crystallography helper functions (:mod:`sknano.core.crystallography._extras`)
==============================================================================

.. currentmodule:: sknano.core.crystallography._extras

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import numpy as np

__all__ = ['supercell_lattice_points']


def supercell_lattice_points(supercell_matrix):
    """Generate the fractional coordinates of lattice points in a \
        supercell lattice.

    Modified implementation of \
    :func:`pymatgen:pymatgen.util.coord_utils.lattice_points_in_supercell`

    Parameters
    ----------
    supercell_matrix : array_like

    Returns
    -------
    :class:`~numpy:numpy.ndarray`
        numpy array of the fractional coordinates of lattice points
    """
    diagonals = np.array(
        [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1],
         [1, 1, 0], [1, 1, 1]])
    d_points = np.dot(diagonals, np.asarray(supercell_matrix))

    mins = np.min(d_points, axis=0)
    maxes = np.max(d_points, axis=0) + 1

    ar = np.arange(mins[0], maxes[0])[:, None] * np.array([1, 0, 0])[None, :]
    br = np.arange(mins[1], maxes[1])[:, None] * np.array([0, 1, 0])[None, :]
    cr = np.arange(mins[2], maxes[2])[:, None] * np.array([0, 0, 1])[None, :]

    all_points = ar[:, None, None] + br[None, :, None] + cr[None, None, :]
    all_points = all_points.reshape((-1, 3))

    frac_points = \
        np.asarray(np.dot(all_points, np.linalg.inv(supercell_matrix)))

    frac_points = frac_points[np.all(frac_points < 1 - 1e-10, axis=1)
                              & np.all(frac_points >= -1e-10, axis=1)]
    assert len(frac_points) == round(abs(np.linalg.det(supercell_matrix)))
    return frac_points
