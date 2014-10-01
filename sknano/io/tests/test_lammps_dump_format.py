# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from pkg_resources import resource_filename
from sknano.io import DUMPData, DUMPReader, DUMPWriter


def test_reader():
    datafile = resource_filename('sknano', 'data/lammpstrj/dump.peptide')
    reader = DUMPReader(datafile)


def test_dumpio():
    datafile = resource_filename('sknano', 'data/lammpstrj/dump.peptide')
    dump = DUMPData(datafile)


if __name__ == '__main__':
    nose.runmodule()
