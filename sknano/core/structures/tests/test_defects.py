# -*- coding: utf-8 -*-
#
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import nose
from nose.tools import *
from sknano.structures import Vacancy, SingleVacancy, DoubleVacancy, \
    TripleVacancy


def test1():
    assert_true(isinstance(SingleVacancy, Vacancy))
    assert_true(isinstance(DoubleVacancy, Vacancy))
    assert_true(isinstance(TripleVacancy, Vacancy))


if __name__ == '__main__':
    nose.runmodule()
