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

from ._base import StructureBase

__all__ = ['Vacancy', 'SingleVacancy', 'DoubleVacancy', 'TripleVacancy',
           'Defect']

Vacancy = namedtuple('Vacancy', ['type', 'size'])
SingleVacancy = Vacancy(type='single', size=1)
DoubleVacancy = Vacancy(type='double', size=2)
TripleVacancy = Vacancy(type='triple', size=3)


class Defect(StructureBase):
    """Base class representation of a structure defect.

    Parameters
    ----------

    """
    def __init__(self):
        pass

    def todict(self):
        """Return :class:`~python:dict` of `Defect` attributes."""
        return dict()
