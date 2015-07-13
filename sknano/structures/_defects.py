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

# import numpy as np

from ._base import StructureBase

__all__ = ['Defect']


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
