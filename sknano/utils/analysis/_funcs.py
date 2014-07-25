# -*- coding: utf-8 -*-
"""
========================================================================
Structure analysis functions (:mod:`sknano.utils.analysis._funcs`)
========================================================================

.. currentmodule:: sknano.utils.analysis._funcs

"""
from __future__ import absolute_import, division, print_function

import numpy as np

__all__ = ['find_target_atom']


def find_target_atom(atoms, target_coords=None, search_radius=2.5,
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
    target_atom :class:`~sknano.structure_io.atoms.Atom`
        An :class:`~sknano.structure_io.atoms.Atom` instance.

    """
    atom_tree = atoms.atom_tree
    print('Starting search radius: {}'.format(search_radius))
    target_atom = None
    while True:
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
            errmsg = 'No target atom found within search radius.\n' + \
                'Increasing search radius by 1 to {}'.format(search_radius)
            print(errmsg)
        else:
            print('Found target atom.')
            return target_atom
