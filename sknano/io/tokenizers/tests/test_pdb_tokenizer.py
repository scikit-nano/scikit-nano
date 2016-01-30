# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import os

import nose
from nose.tools import assert_equal
# from pkg_resources import resource_filename
from sknano.io import PDBReader  # , PDBData, PDBWriter, PDB2DATAConverter
from sknano.testing import IOTestFixture


class PDBTestFixture(IOTestFixture):
    pass


class Tests(PDBTestFixture):

    pass


if __name__ == '__main__':
    nose.runmodule()
