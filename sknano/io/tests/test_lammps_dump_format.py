# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from pkg_resources import resource_filename
from sknano.io import DUMPData, DUMPReader, DUMPWriter


def test_reader():
    datafile = resource_filename('sknano', 'data/nanotubes/1010_1cell.data')
    reader = DUMPReader(datafile)
    atoms = reader.atoms
    assert_not_equal(atoms.Natoms, 20)
    assert_equal(atoms.Natoms, 40)


if __name__ == '__main__':
    nose.runmodule()
