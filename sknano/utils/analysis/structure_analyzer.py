# -*- coding: utf-8 -*-
"""
===============================================================================
Structure analyzer class (:mod:`sknano.utils.analysis.structure_analyzer`)
===============================================================================

.. currentmodule:: sknano.utils.analysis.structure_analyzer

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

__all__ = ['StructureAnalyzer']


class StructureAnalyzer:
    """Class for structure analysis.

    Parameters
    ----------
    structure_data : str
        structure data file

    """
    def __init__(self, structure, **kwargs):
        self.structure = structure
        self.structure.atoms.update_attrs()

    def analyze_POAV(self):
        self.structure.atoms.compute_POAVs()
