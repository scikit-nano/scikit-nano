# -*- coding: utf-8 -*-
"""
===============================================================================
Core package code (:mod:`sknano.core`)
===============================================================================

.. currentmodule:: sknano.core

Contents
========
Core code for package development and general use.

Data structures & algorithms
-----------------------------
.. autosummary::
   :toctree: generated/

   dedupe
   rezero_array


Iterator functions
------------------

.. autosummary::
   :toctree: generated/

   cyclic_pairs
   take
   tabulate
   consume
   nth
   quantify
   padnone
   ncycles
   dotproduct
   flatten
   repeatfunc
   pairwise
   grouper
   roundrobin
   partition
   powerset
   unique_elements
   iter_except
   first_true
   random_product
   random_permutation
   random_combination
   random_combination_with_replacement


Meta functions/classes
-----------------------

.. autosummary::
   :toctree: generated/

   check_type
   deprecated
   get_object_signature
   memoize
   timethis
   method_func
   with_doc
   make_sig
   ClassSignature


I/O functions
--------------

.. autosummary::
   :toctree: generated/

   get_fname
   get_fpath

String functions
-----------------

.. autosummary::
   :toctree: generated/

   pluralize

Custom container datatypes
---------------------------

.. autosummary::
   :toctree: generated/

   ListBasedSet
   UserList
   frozendict

Sub-packages
-------------

* atoms (:mod:`sknano.core.atoms`)
* crystallography (:mod:`sknano.core.crystallography`)
* geometric regions (:mod:`sknano.core.geometric_regions`)
* math (:mod:`sknano.core.math`)
* molecules (:mod:`sknano.core.molecules`)
* physics (:mod:`sknano.core.physics`)
* refdata (:mod:`sknano.core.refdata`)

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from ._collections import *
from ._extras import *
from ._io import *
from ._itertools import *
from ._meta import *
from ._strings import *

__all__ = [s for s in dir() if not s.startswith('_')]
