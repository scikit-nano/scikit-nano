# -*- coding: utf-8 -*-
"""
========================================================================
Helper functions for structure analysis (:mod:`sknano.analysis._funcs`)
========================================================================

.. currentmodule:: sknano.analysis._funcs

"""
from __future__ import absolute_import, division, print_function

import numpy as np

__all__ = ['find_target_atom']


def find_target_atom(atoms, target_coords=None, search_radius=float,
                     nearest_target=False):
    """Search for atom closest to target location.

    Parameters
    ----------
    atoms : :class:`~sknano.structure_io.atoms.Atoms`
        An :class:`~sknano.structure_io.atoms.Atoms` instance.
    target_coords : array_like
        An array or list of :math:`x,y,z` coordinates
    search_radius : float
        Cutoff distance to search within in units of
        :class:`~sknano.structure_io.atoms.Atoms` style.
    nearest_target : bool, optional

    Returns
    -------
    :class:`~sknano.structure_io.atoms.Atom`
        An :class:`~sknano.structure_io.atoms.Atom` instance.

    """
    atom_tree = atoms.atom_tree
    target_index = None

    if nearest_target:
        target_index = \
            atom_tree.query(target_coords,
                            distance_upper_bound=search_radius)[-1]
    else:
        target_index = \
            np.random.choice(
                atom_tree.query_ball_point(target_coords, search_radius))

    return atoms[target_index]
