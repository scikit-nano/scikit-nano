# -*- coding: utf-8 -*-
"""
===============================================================================
Core code for package development and general use (:mod:`sknano.core`)
===============================================================================

.. currentmodule:: sknano.core

Contents
========

Functions
----------

Custom iterators:

.. autosummary::
   :toctree: generated/

   cyclic_pairs

Array functions:

.. autosummary::
   :toctree: generated/

   rezero_array

Meta functions:

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

I/O functions:

.. autosummary::
   :toctree: generated/

   get_fname
   get_fpath

String functions:

.. autosummary::
   :toctree: generated/

   pluralize

Classes
-------

.. autosummary::
   :toctree: generated/

   UserList

Sub-packages
-------------

* atoms (:mod:`sknano.core.atoms`)
* math (:mod:`sknano.core.math`)
* physics (:mod:`sknano.core.physics`)
* refdata (:mod:`sknano.core.refdata`)

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from ._extras import *
from ._io import *
from ._meta import *
from ._strings import *
from ._user_list import *

__all__ = [s for s in dir() if not s.startswith('_')]
