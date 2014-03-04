# -*- coding: utf-8 -*-
"""
=======================================================================
Helper functions, LUTs, reference data, etc. (:mod:`sknano.tools`)
=======================================================================

.. currentmodule:: sknano.tools

Contents
========

Abstract data structures
-------------------------
.. autosummary::
   :toctree: generated/

   Point
   Vector
   Quaternion

Helper functions for chirality data
------------------------------------
.. autosummary::
   :toctree: generated/

   cmp_Ch
   filter_Ch
   filter_Ch_list
   generate_Ch_list
   get_Ch_indices
   get_Ch_type

'Core' functions that don't fit under any specific category
------------------------------------------------------------
.. autosummary::
   :toctree: generated/

   check_type

Helper functions for I/O
-------------------------
.. autosummary::
   :toctree: generated/

   get_fpath

Helper functions for manipulating strings
------------------------------------------
.. autosummary::
   :toctree: generated/

   plural_word_check

Math functions
---------------
.. autosummary::
   :toctree: generated/

   totient_func

Helper functions for linear algebra transforms
-----------------------------------------------
.. autosummary::
   :toctree: generated/

   rotation_matrix

LUTs (Look-up 'tables' - lists and dictionaries)
-------------------------------------------------
.. autosummary::
   :toctree: generated/

   chiral_type_name_mappings
   xyz_axes

Sub-packages
============
.. autosummary::
   :toctree: generated/

   refdata

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from ._chirality_funcs import *
from ._corefuncs import *
from ._coremath import *
from ._iofuncs import *
from ._luts import *
from ._mathfuncs import *
from ._strfuncs import *
from ._transforms import *

__all__ = [s for s in dir() if not s.startswith('_')]
