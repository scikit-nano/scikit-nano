# -*- coding: utf-8 -*-
"""
==============================================================================
Defect structure classes (:mod:`sknano.structures._defects`)
==============================================================================

.. currentmodule:: sknano.structures._defects

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from collections import namedtuple

# import numpy as np

# from sknano.core import frozendict

# from ._base import NanoStructureBase

__all__ = ['Vacancy', 'SingleVacancy', 'DoubleVacancy', 'TripleVacancy',
           'PointDefect']

Vacancy = namedtuple('Vacancy', ['type', 'size'])
SingleVacancy = Vacancy(type='single', size=1)
DoubleVacancy = Vacancy(type='double', size=2)
TripleVacancy = Vacancy(type='triple', size=3)

# VACANCY_TYPE = frozendict({'single': 1, 'double': 2, 'triple': 3})
# __all__ += ['VACANCY_TYPE']


class PointDefect:
    """Base class representation of point defects.

    Parameters
    ----------

    """
    def __init__(self):
        pass

    def todict(self):
        """Return :class:`~python:dict` of `Defect` attributes."""
        return dict()


class VacancyDefect(PointDefect):
    pass


class FrenkelDefect(PointDefect):
    pass


class Impurities(PointDefect):
    pass


class InterstitialDefect(PointDefect):
    pass


class TopologicalDefect:
    """Base class representation of topological defects."""
    pass


class StoneWalesDefect(TopologicalDefect):
    pass


class Adatom(TopologicalDefect):
    pass


class ComplexDefect:
    pass
