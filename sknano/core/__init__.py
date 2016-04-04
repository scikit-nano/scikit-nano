# -*- coding: utf-8 -*-
"""
===============================================================================
Core package code (:mod:`sknano.core`)
===============================================================================

.. currentmodule:: sknano.core

Contents
========
Core code for package development and general use.

Sub-packages
-------------

* analysis (:mod:`sknano.core.analysis`)
* atoms (:mod:`sknano.core.atoms`)
* crystallography (:mod:`sknano.core.crystallography`)
* geometric regions (:mod:`sknano.core.geometric_regions`)
* math (:mod:`sknano.core.math`)
* physics (:mod:`sknano.core.physics`)
* refdata (:mod:`sknano.core.refdata`)
* structures (:mod:`sknano.core.structures`)

Data structures & algorithms
-----------------------------
.. autosummary::
   :toctree: generated/

   dedupe
   rezero_array
   minmax


Iterator functions
------------------
.. autosummary::
   :toctree: generated/

   cyclic_pairs
   take
   tabulate
   tail
   consume
   nth
   all_equal
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
   unique_everseen
   unique_justseen
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
   deprecate_kwarg
   get_object_signature
   lazy_property
   logged
   memoize
   optional_debug
   timethis
   typeassert
   typed_property
   make_sig
   BaseClass
   Cached
   ClassSignature
   NoInstances
   Singleton
   method_func
   with_doc


I/O functions
--------------
.. autosummary::
   :toctree: generated/

   get_fname
   get_fpath
   listdir_dirnames
   listdir_fnames
   listdir
   loadobj
   dumpobj

String functions/classes
-------------------------
.. autosummary::
   :toctree: generated/

   asbool
   asdict
   aslist
   asset
   astuple
   asint
   asfloat
   map_function
   map_operator
   ordinal_form
   pluralize
   TabulateMixin
   obj_mro_str

Custom container datatypes
---------------------------

.. autosummary::
   :toctree: generated/

   ListBasedSet
   UserList
   frozendict

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
