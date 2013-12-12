# -*- coding: utf-8 -*-
"""
=========================================================================
Structure data (:mod:`sknano.structure_io._structure_data`)
=========================================================================

.. currentmodule:: sknano.structure_io._structure_data

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from collections import OrderedDict

__all__ = ['StructureData', 'LAMMPSDATA']


class StructureData(object):
    """Base class defining common properties for structure formats."""

    def __init__(self):
        self._properties = OrderedDict()

    def properties(self):
        """OrderedDict of format properties."""
        return self._properties


class LAMMPSDATA(StructureData):
    """
    Class defining the structure file format for LAMMPS data.

    """
    def __init__(self):
        super(LAMMPSDATA, self).__init__()
        pass
