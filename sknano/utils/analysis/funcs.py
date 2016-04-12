# -*- coding: utf-8 -*-
"""
========================================================================
Structure analysis functions (:mod:`sknano.utils.analysis.funcs`)
========================================================================

.. currentmodule:: sknano.utils.analysis.funcs

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from collections import deque

import numpy as np

from sknano.core.atoms import Atoms


__all__ = ['find_defect_chains', 'find_target_atom']


def find_defect_chains(atoms, defect_condition, max_length=None,
                       exclude=None, include_root=True, **kwargs):
    """Analyze atoms for chains of defects."""
    if not isinstance(atoms, Atoms):
        raise TypeError('Expected an `Atoms` object')

    if max_length is None:
        max_length = len(atoms)
    if exclude is None:
        exclude = []

    chains = []
    seen = []
    for atom in atoms:
        chain = []
        queue = deque([natom for natom in atom.neighbors])
        while queue and len(chain) <= max_length:
            natom = queue.popleft()
            if defect_condition(natom) and \
                    all([natom not in atoms_ for atoms_ in
                         (atoms, chain, chains, seen, exclude)]):
                chain.append(natom)
                queue.extend([nnatom for nnatom in natom.neighbors if
                              all([nnatom not in atoms_ for atoms_ in
                                   (atoms, chain, chains, seen, exclude)])])
            seen.append(natom)
        if len(chain) > 0:
            if include_root:
                chain.insert(0, atom)
            defect_atoms = \
                atoms.__class__(chain, update_item_class=False,
                                **atoms.kwargs)
            chains.append(defect_atoms)

    return chains


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
    target_atom : :class:`~sknano.core.atoms.Atom`
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
