# -*- coding: utf-8 -*-
"""
======================================================================
Meta functions and class (:mod:`sknano.core._meta`)
======================================================================

.. currentmodule:: sknano.core._meta

Contents
========

"""
from __future__ import absolute_import, division, print_function

from functools import wraps
import inspect
import sys
import time
import warnings
#warnings.resetwarnings()

import numpy as np
from numpy.compat import formatargspec, getargspec

__all__ = ['check_type', 'deprecated', 'get_object_signature', 'memoize',
           'method_function', 'methodfunc', 'removed_package_warning',
           'timethis', 'with_doc']


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


def deprecated(replacement=None):
    """Deprecated decorator.

    A decorator which can be used to mark functions as deprecated.
    replacement is a callable that will be called with the same args
    as the decorated function.

    Source: http://code.activestate.com/recipes/577819-deprecated-decorator/
    Original Author: Giampaolo Rodola' <g.rodola [AT] gmail [DOT] com>
    License: MIT

    Parameters
    ----------
    replacement : callable, optional

    Returns
    -------
    decorator : decorator function

    """
    def decorate(func):
        msg = "{!s} is deprecated".format(func.__name__)
        if replacement is not None:
            msg += "; use {!s} instead".format(replacement)
        if func.__doc__ is None:
            func.__doc__ = msg

        @wraps(func)
        def wrapper(*args, **kwargs):
            warnings.warn_explicit(
                "Call to deprecated function {}".format(func.__name__),
                category=DeprecationWarning,
                filename=func.__code__.co_filename,
                lineno=func.__code__.co_firstlineno + 1
            )
            return func(*args, **kwargs)
        return wrapper
    return decorate


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
        key = (f, tuple(args), frozenset(list(kwargs.items())))
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


def removed_package_warning(oldpkg, newpkg=None):
    msg = '\n\n#TODO: UPDATE THIS WARNING MESSAGE.\n'
    #msg = '\n\n{:<74.74}\n'.format('{:*^80}'.format(' NB '))
    #    "As of version 0.3, the {!r} package was ".format(oldpkg)
    #if newpkg is None:
    #    msg += "removed.\n"
    #else:
    #    msg += "renamed.\n\n" + \
    #        "Replace imports from: {!r}\n".format(oldpkg) + \
    #        "with imports from: {!r}\n".format(newpkg)
    #msg += "Also review the module specific Warnings issued by\n" + \
    #    "module imports from a deprecated package location, as\n" + \
    #    "in some cases the module code may have moved to a different\n" + \
    #    "package entirely.\n\n"
    #msg += "Please update your code accordingly as these warnings and\n" + \
    #    "associated fallback code will likely be removed in a future version."
    #msg += "\n{}\n".format(74 * '*')

    warnings.warn(msg, ImportWarning)


def timethis(func):
    """Decorator that reports execution time."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(func.__name__, end-start)
        return result
    return wrapper


class with_doc(object):
    """Decorator class to combine class method docstrings.

    Parameters
    ----------
    method : method object

    """
    def __init__(self, method, use_header=True):
        self.method = method
        if use_header:
            self.header = "\n{}Notes\n{}\n".format(4 * ' ', 4 * ' ', 4 * ' ')
        else:
            self.header = ''

    def __call__(self, new_method):
        new_doc = new_method.__doc__
        original_doc = self.method.__doc__
        header = self.header

        if original_doc and new_doc:
            new_method.__doc__ = \
                '\n'.join([''.join([4 * ' ', s]) for s in
                           (original_doc, header, new_doc)])
        elif original_doc:
            new_method.__doc__ = original_doc

        return new_method
