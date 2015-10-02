#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *

# import numpy as np

from sknano.utils.analysis import find_target_atom
from sknano.testing import generate_atoms


def test_find_target_atom():
    atoms = generate_atoms(generator_class='SWNTGenerator', n=10, m=10, nz=3)
    atoms.center_centroid()
    atoms.assign_unique_ids()
    atoms.update_attrs()
    target_atom = find_target_atom(atoms, target_coords=[0.0, 5.0, 0.0],
                                   search_radius=1.0, nearest_target=True)

    print('atoms.z.min(), atoms.z.max(): {}, {}'.format(
        atoms.z.min(), atoms.z.max()))
    print('atom_ids: {}'.format(atoms.atom_ids))
    print('CNs: {}'.format(atoms.coordination_numbers))

    print('target_atom: {}'.format(target_atom))


if __name__ == '__main__':
    nose.runmodule()
