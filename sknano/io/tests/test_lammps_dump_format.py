# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function

import nose
from nose.tools import *
from pkg_resources import resource_filename
from sknano.io import DUMPData, DUMPReader, DUMPWriter


#def test_reader():
#    dumpfile = resource_filename('sknano', 'data/lammpstrj/dump.peptide')
#    dumpdata = DUMPReader(dumpfile)
#    print('dumpdata.nsnaps: {}'.format(dumpdata.nsnaps))


def test_dump():
    #dumpfile = resource_filename('sknano', 'data/lammpstrj/dump.peptide')
    dumpfile = resource_filename('sknano', 'data/lammpstrj/1010+ion.dump')
    dumpdata = DUMPData(dumpfile)
    #print('dumpdata.nsnaps: {}'.format(dumpdata.nsnaps))


if __name__ == '__main__':
    nose.runmodule()
