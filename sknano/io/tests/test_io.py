# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from pprint import pprint

import unittest

from sknano.io import DATAData, DATAReader, DATAWriter, DATA2XYZConverter, \
    XYZData, XYZReader, XYZWriter, XYZ2DATAConverter


class TestDATAIO(unittest.TestCase):

    def setUp(self):
        self.infile = '1010r_1cell.data'

    def test_reader(self):
        reader = DATAReader(self.infile)
        atoms = reader.atoms
        print('Natoms: {}'.format(atoms.Natoms))
        print('atoms: {}'.format(atoms))


class TestXYZIO(unittest.TestCase):

    def setUp(self):
        self.infile = '1010r_1cell.xyz'

    def test_reader(self):
        reader = XYZReader(self.infile)
        atoms = reader.atoms
        print('Natoms: {}'.format(atoms.Natoms))
        print('atoms: {}'.format(atoms))


if __name__ == '__main__':
    unittest.main()
