# -*- coding: utf-8 -*-
"""
=======================================================================
Tools for analysis, helper functions, LUTs, etc. (:mod:`sknano.tools`)
=======================================================================

.. currentmodule:: sknano.tools

.. versionadded:: 0.2.6
   `GrapheneVacancyGenerator` and `NanotubeVacancyGenerator` implemented

Contents
========

Tools for manipulating nano-structures
--------------------------------------

.. autosummary::
   :toctree: generated/

   VacancyGenerator
   GrapheneVacancyGenerator
   NanotubeVacancyGenerator

Helper Functions
-----------------

.. autosummary::
   :toctree: generated/

   cmp_Ch
   generate_Ch_list
   get_Ch_indices
   get_Ch_type

LUTs (Look-up 'tables' - lists and dictionaries)
-------------------------------------------------

.. autosummary::
   :toctree: generated/

   chiral_type_name_mappings
   xyz_axes

Custom exception classes for handling errors
--------------------------------------------

.. autosummary::
   :toctree: generated/

   VacancyGeneratorError

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from ._funcs import *
from ._luts import *
from ._vacancy_generator import *

__all__ = [s for s in dir() if not s.startswith('_')]
