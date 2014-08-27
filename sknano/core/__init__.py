# -*- coding: utf-8 -*-
"""
===============================================================================
Core code for testing, development, and general use (:mod:`sknano.core`)
===============================================================================

.. currentmodule:: sknano.core

Contents
========

Core package helper functions
------------------------------
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

Functions for file I/O
-----------------------
.. autosummary::
   :toctree: generated/

   get_fpath

Functions for string manipulation
----------------------------------
.. autosummary::
   :toctree: generated/

   plural_word_check

Sub-packages
------------
.. autosummary::
   :toctree: generated/

   atoms
   math
   refdata
   testing

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._core import *
from ._decorators import *
from ._iofuncs import *
from ._luts import *
from ._strfuncs import *
from ._warnings import *

__all__ = [s for s in dir() if not s.startswith('_')]
