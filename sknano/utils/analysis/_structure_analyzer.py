# -*- coding: utf-8 -*-
"""
===============================================================================
Structure analyzer class (:mod:`sknano.utils.analysis._structure_analyzer`)
===============================================================================

.. currentmodule:: sknano.utils.analysis._structure_analyzer

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

__all__ = ['StructureAnalyzer']

#import numpy as np


class StructureAnalyzer:
    """Class for structure analysis.

    Parameters
    ----------
    structure_data : str
        structure data file

    """
    def __init__(self, structure_data, **kwargs):
        self.structure_data = structure_data
        self.structure_data.atoms.update_attrs()

    def analyze_POAV(self):
        self.structure_data.atoms.compute_POAVs()
