# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import os

import nose
from nose.tools import assert_equal, assert_raises
# from pkg_resources import resource_filename
from sknano.io import PDBReader, StructureIOError
from sknano.testing import IOTestFixture


class PDBTestFixture(IOTestFixture):
    def setUp(self):
        atoms = self.atoms = self.pdb_reader.atoms
        atoms.assign_unique_ids()
        atoms.assign_unique_types()
        atoms.update_attrs()

    @property
    def pdbfile1(self):
        return os.path.join(os.path.dirname(__file__), 'data/pdb1iar.ent')

    @property
    def pdbfile2(self):
        return os.path.join(os.path.dirname(__file__), 'data/Nanotube.pdb')


class Tests(PDBTestFixture):

    def test1(self):
        atoms = self.atoms
        assert_equal(atoms.Natoms, 40)
        assert_equal(list(set(atoms.atom_ids)), list(range(1, 41)))

    def test2(self):
        atoms = PDBReader(self.pdbfile1).atoms
        # print('atoms.Natoms: {}'.format(atoms.Natoms))
        assert_equal(2795, atoms.Natoms)

    def test3(self):
        with assert_raises(StructureIOError):
            PDBReader(self.pdbfile2)


if __name__ == '__main__':
    nose.runmodule()
