# -*- coding: utf-8 -*-
"""
========================================================================
Structure analysis functions (:mod:`sknano.utils.analysis._funcs`)
========================================================================

.. currentmodule:: sknano.utils.analysis._funcs

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import numpy as np

__all__ = ['find_target_atom']


def find_target_atom(atoms, target_coords=None, search_radius=1.0,
                     nearest_target=False, max_search_radius=None):
    """Search for atom closest to target location.

    Parameters
    ----------
    atoms : :class:`~sknano.core.atoms.Atoms`
        An :class:`~sknano.core.atoms.Atoms` instance.
    target_coords : array_like
        An array or list of :math:`x,y,z` coordinates
    search_radius : :class:`~python:float`, optional
        Search radius
    nearest_target : :class:`~python:bool`, optional
    max_search_radius : :class:`~python:float`, optional

    Returns
    -------
    target_atom :class:`~sknano.core.atoms.Atom`
        An :class:`~sknano.core.atoms.Atom` instance.

    """
    print('Starting search radius: {}'.format(search_radius))
    target_atom = None
    if max_search_radius is None:
        max_search_radius = 100. * search_radius
    atom_tree = atoms.atom_tree
    while search_radius <= max_search_radius:
        try:
            target_index = None
            if nearest_target:
                target_index = \
                    atom_tree.query(target_coords,
                                    distance_upper_bound=search_radius)[-1]
            else:
                target_index = \
                    np.random.choice(
                        atom_tree.query_ball_point(target_coords,
                                                   search_radius))
            target_atom = atoms[target_index]
        except (IndexError, ValueError):
            search_radius += 1
            print('No target atom within search radius.\n'
                  'Increasing radius by 1 unit to {}'.format(search_radius))
        else:
            print('Found target atom.')
            return target_atom
    else:
        print('Maximum search radius exceeded. No target atom found.')
        return None
