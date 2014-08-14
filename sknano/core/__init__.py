# -*- coding: utf-8 -*-
"""
======================================================================
Core package modules (:mod:`sknano.core`)
======================================================================

.. currentmodule:: sknano.core

Contents
========

Abstract data structures for Atom objects
-----------------------------------------
.. autosummary::
   :toctree: generated/

   Atom
   Atoms

Core helper functions
----------------------
.. autosummary::
   :toctree: generated/

   check_type
   get_object_signature
   memoize
   method_function
   methodfunc

Decorator functions
--------------------
.. autosummary::
   :toctree: generated/

   deprecated
   with_doc

Helper functions for I/O
-------------------------
.. autosummary::
   :toctree: generated/

   get_fpath

Math functions
---------------
.. autosummary::
   :toctree: generated/

   totient_func

Abstract mathematical data structures
--------------------------------------
.. autosummary::
   :toctree: generated/

   Point
   Vector

Helper functions for manipulating strings
------------------------------------------
.. autosummary::
   :toctree: generated/

   plural_word_check

Helper functions for linear algebra transforms
-----------------------------------------------
.. autosummary::
   :toctree: generated/

   rotate_point
   rotation_matrix
   transformation_matrix

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._atom import *
from ._atoms import *
from ._core import *
from ._decorators import *
from ._iofuncs import *
from ._luts import *
from ._mathfuncs import *
from ._npcoremath import *
from ._strfuncs import *
from ._transforms import *
from ._warnings import *

__all__ = [s for s in dir() if not s.startswith('_')]
