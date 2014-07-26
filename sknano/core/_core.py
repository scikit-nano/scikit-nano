# -*- coding: utf-8 -*-
"""
========================================================
Core classes and functions (:mod:`sknano.core._core`)
========================================================

.. currentmodule:: sknano.core._core

"""
from __future__ import absolute_import, division, print_function

from functools import wraps
import inspect
import sys

import numpy as np
from numpy.compat import formatargspec, getargspec

__all__ = ['check_type', 'get_object_signature', 'memoize',
           'method_function', 'methodfunc']


def check_type(value, allowed_types=()):
    """Check object `value` type against tuple of `allowed_types`.

    Parameters
    ----------
    value : `object`
    allowed_types : tuple
        tuple of allowed classes and/or types

    Raises
    ------
    `TypeError`
        If `value` fails `isinstance` check against `allowed_types`.

    """
    if not isinstance(value, allowed_types):
        raise TypeError('{} does not have an allowed type.\n'.format(value) +
                        '(Allowed Type(s): {})'.format(allowed_types))


def get_object_signature(obj):
    """
    Get the signature from obj
    """
    try:
        sig = formatargspec(*getargspec(obj))
    except TypeError:
        sig = ''

    return sig


def memoize(f, cache={}):
    """memoization function to cache dict mapping"""
    @wraps(f)
    def g(*args, **kwargs):
        key = (f, tuple(args), frozenset(kwargs.items()))
        if key not in cache:
            cache[key] = f(*args, **kwargs)
        return cache[key].copy()
    return g


def method_function(module, classobj):
    # get the module as an object
    moduleobj = sys.modules[module]

    # Iterate over the methods of the class and dynamically create a function
    # for each method that calls the method and add it to the current module
    for member in inspect.getmembers(classobj, predicate=inspect.ismethod):
        method = member[0]

        # get the bound method
        func = getattr(classobj, method)

        # add the function to the current module
        setattr(moduleobj, method, func)


class methodfunc(object):
    """Define functions from existing class methods.

    This class is based off of the the `numpy`
    :class:`~numpy:numpy.ma._frommethod` class.

    Parameters
    ----------
    classname : str
        Name of the class containing `methodname`
    methodname : str
        Name of the method to transform.

    """
    def __init__(self, classobj, methodname, reversed=False):
        self._classobj = classobj
        self.__name__ = methodname
        self.__doc__ = self.getdoc()
        self.reversed = reversed

    def getdoc(self):
        "Return the doc of the function (from the doc of the method)."
        meth = getattr(self._classobj, self.__name__, None) or \
            getattr(np, self.__name__, None)
        signature = self.__name__ + get_object_signature(meth)
        if meth is not None:
            doc = """    {}\n{}""".format(signature,
                                          getattr(meth, '__doc__', None))
            return doc

    def __call__(self, a, *args, **params):
        if self.reversed:
            args = list(args)
            arr = args[0]
            args[0] = a
            a = arr
        # Get the method from the array (if possible)
        method_name = self.__name__
        method = getattr(a, method_name, None)
        if method is not None:
            return method(*args, **params)
        # Still here ? Then a is not a classobj
        method = getattr(self._classobj, method_name, None)
        if method is not None:
            return method(self._classobj(a), *args, **params)
        # Still here ? OK, let's call the corresponding np function
        method = getattr(np, method_name)
        return method(a, *args, **params)
