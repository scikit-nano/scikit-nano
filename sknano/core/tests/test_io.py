# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename
import nose
from nose.tools import *

from sknano.core import get_fname, get_fpath, listdir_dirnames, \
    listdir_fnames, listdir


def test1():
    datadir = resource_filename('sknano', 'data')
    dirnames, fnames = listdir(datadir)
    assert_true('__init__.py' in fnames)
    assert_true('fullerenes' in dirnames)
    assert_equal(dirnames, listdir_dirnames(datadir))
    assert_equal(fnames, listdir_fnames(datadir))


def test2():
    datadir = resource_filename('sknano', 'data')
    filterfunc = lambda name: '_' not in name
    dirnames, fnames = listdir(datadir, filterfunc=filterfunc,
                               filter_fnames=True, filter_dirnames=True)
    assert_true('__pycache__' not in dirnames)
    assert_true('__init__.py' not in fnames)
    assert_equal(dirnames,
                 listdir_dirnames(datadir, filterfunc=filterfunc))
    assert_equal(fnames,
                 listdir_fnames(datadir, filterfunc=filterfunc))


if __name__ == '__main__':
    nose.runmodule()
