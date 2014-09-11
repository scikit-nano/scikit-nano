# -*- coding: utf-8 -*-
"""
===============================================================================
Core package code for development and general use (:mod:`sknano.core`)
===============================================================================

.. currentmodule:: sknano.core

Contents
========

Core meta functions/classes
------------------------------
.. autosummary::
   :toctree: generated/

   check_type
   deprecated
   get_object_signature
   memoize
   method_function
   timethis

   methodfunc
   with_doc

I/O functions
--------------
.. autosummary::
   :toctree: generated/

   get_fname
   get_fpath

Functions for string manipulation
----------------------------------
.. autosummary::
   :toctree: generated/

   pluralize

Sub-packages
------------
.. autosummary::
   :toctree: generated/

   atoms
   math
   refdata

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._extras import *
from ._io import *
from ._meta import *
from ._strings import *

__all__ = [s for s in dir() if not s.startswith('_')]
