# -*- coding: utf-8 -*-
"""
==============================================================================
Structure analysis module (:mod:`sknano.utils.analysis._structure_analyzer`)
==============================================================================

.. currentmodule:: sknano.utils.analysis._structure_analyzer

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

__all__ = ['StructureAnalyzer']


class StructureAnalyzer(object):
    u"""Class for structure analysis.

    Parameters
    ----------
    structure_data : str
        structure data file

    """
    def __init__(self, structure_data, **kwargs):
        self.structure_data = structure_data
