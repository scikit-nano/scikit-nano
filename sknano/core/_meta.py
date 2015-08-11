# -*- coding: utf-8 -*-
"""
======================================================================
Meta functions and classes (:mod:`sknano.core._meta`)
======================================================================

.. currentmodule:: sknano.core._meta

Contents
========

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from abc import ABCMeta, abstractmethod
from inspect import signature, Signature, Parameter
from functools import wraps, partial
import inspect
import logging
import time
import warnings
import weakref
# warnings.resetwarnings()

import numpy as np
from numpy.compat import formatargspec, getargspec

__all__ = ['attach_wrapper', 'logged', 'check_type', 'deprecated',
           'get_object_signature', 'memoize',
           'method_func', 'optional_debug', 'removed_package_warning',
           'timethis', 'typeassert', 'typed_property', 'with_doc',
           'make_sig', 'ClassSignature', 'Cached', 'NoInstances',
           'Singleton', 'BaseClass']


def attach_wrapper(obj, func=None):
    """Utility decorator to attach a function as an attribute of `obj`.

    Parameters
    ----------
    obj : :class:`~python:object`
    func : {None, callable}

    """
    if func is None:
        return partial(attach_wrapper, obj)
    setattr(obj, func.__name__, func)
    return func


def logged(func=None, *, level=logging.DEBUG, name=None, message=None):
    """Decorator to add logging to a function.

    Parameters
    ----------
    func : callable
        Decorated function
    level : :class:`~python:int`
        Logging level
    name : :class:`~python:str`, optional
        Logger name
    message : :class:`~python:str`, optional
        Log message
    """
    if func is None:
        return partial(logged, level=level, name=name, message=message)

    logname = name if name else func.__module__
    log = logging.getLogger(logname)
    logmsg = message if message else func.__name__

    @wraps(func)
    def wrapper(*args, **kwargs):
        log.log(level, logmsg)
        return func(*args, **kwargs)

    # Attach setter functions
    @attach_wrapper(wrapper)
    def set_level(newlevel):
        nonlocal level
        level = newlevel

    @attach_wrapper(wrapper)
    def set_message(newmsg):
        nonlocal logmsg
        logmsg = newmsg

    return wrapper


def check_type(obj, allowed_types=()):
    """Check object type against tuple of `allowed_types`.

    Parameters
    ----------
    obj : `object`
    allowed_types : tuple
        tuple of allowed classes and/or types

    Raises
    ------
    `TypeError`
        If `value` fails `isinstance` check against `allowed_types`.

    """
    if not isinstance(obj, allowed_types):
        raise TypeError('{} does not have an allowed type.\n'.format(obj) +
                        '(Allowed type(s): {})'.format(allowed_types))


def deprecated(replacement=None):
    """Decorator to mark functions as deprecated.

    Parameters
    ----------
    replacement : callable, optional
        callable that will be called with the same args
        as the decorated function.

    Returns
    -------
    decorate : decorated function

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


def optional_debug(func):
    """Decorator that adds an optional `debug` keyword argument."""
    if 'debug' in inspect.getargspec(func).args:
        raise TypeError('debug argument already defined')

    @wraps(func)
    def wrapper(*args, debug=False, **kwargs):
        if debug:
            print('Calling', func.__name__)
        return func(*args, **kwargs)

    sig = inspect.signature(func)
    params = list(sig.parameters.values())
    params.append(inspect.Parameter('debug',
                                    inspect.Parameter.KEYWORD_ONLY,
                                    default=False))
    wrapper.__signature__ = sig.replace(parameters=params)
    return wrapper


def removed_package_warning(oldpkg, newpkg=None):
    msg = '\n\n#TODO: UPDATE THIS WARNING MESSAGE.\n'
    # msg = '\n\n{:<74.74}\n'.format('{:*^80}'.format(' NB '))
    #    "As of version 0.3, the {!r} package was ".format(oldpkg)
    # if newpkg is None:
    #    msg += "removed.\n"
    # else:
    #    msg += "renamed.\n\n" + \
    #        "Replace imports from: {!r}\n".format(oldpkg) + \
    #        "with imports from: {!r}\n".format(newpkg)
    # msg += "Also review the module specific Warnings issued by\n" + \
    #    "module imports from a deprecated package location, as\n" + \
    #    "in some cases the module code may have moved to a different\n" + \
    #    "package entirely.\n\n"
    # msg += "Please update your code accordingly as these warnings and\n" + \
    #    "associated fallback code will likely be removed in a future version."
    # msg += "\n{}\n".format(74 * '*')

    warnings.warn(msg, ImportWarning)


def timethis(func):
    """Decorator that reports execution time."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(func.__name__, end - start)
        return result
    return wrapper


def typeassert(*type_args, **type_kwargs):
    """Decorator that enforces type checking of function arguments."""
    def decorate(func):
        # If in optimized mode, disable type checking
        if not __debug__:
            return func

        # Map function argument names to supplied types
        sig = signature(func)
        bound_types = sig.bind_partial(*type_args, **type_kwargs).arguments

        @wraps(func)
        def wrapper(*args, **kwargs):
            bound_values = sig.bind(*args, **kwargs)
            # Enforce type assertions across supplied arguments
            for name, value in bound_values.arguments.items():
                if name in bound_types:
                    if not isinstance(value, bound_types[name]):
                        raise TypeError('Argument {} must be {}'.format(
                                        name, bound_types[name]))
            return func(*args, **kwargs)
        return wrapper
    return decorate


def typed_property(name, expected_type):
    _name = '_' + name

    @property
    def prop(self):
        return getattr(self, _name)

    @prop.setter
    def prop(self, value):
        check_type(value, expected_type)
        setattr(self, _name, value)

    return prop


class BaseClass(metaclass=ABCMeta):
    """ABC defining a common set of attributes/methods for other base classes.

    Attributes
    ----------
    fmtstr

    """

    def __init__(self, *args, verbose=False, debug=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.verbose = verbose
        self.debug = debug
        self.fmtstr = ""

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.fmtstr.format(**self.todict()))

    @property
    def fmtstr(self):
        """Format string.

        User defined `format string`_ that should be set by all
        :class:`BaseClass` sub-classes. The format string should
        contain replacement fields that are the named keyword arguments
        contained in the :class:`~python:dict` returned by the sub-class
        implementation of the :meth:`~BaseClass.todict` method,
        which is required to be overridden by any callable sub-class
        of `BaseClass`.

        .. _Format string: https://docs.python.org/3/library/string.html#format-string-syntax
        """
        return self._fmtstr

    @fmtstr.setter
    def fmtstr(self, value):
        self._fmtstr = value

    @abstractmethod
    def todict(self):
        """Return :class:`~python:dict` of constructor parameters.

        The :class:`~python:dict` should contain the same named
        keyword arguments defined in the replacement fields of the
        :attr:`~BaseClass.fmtstr` defined by any subclass of `BaseClass`.
        """
        return NotImplementedError


class Cached(type):
    """Cached class type."""
    def __init__(self, *args, **kwargs):
        # print('in Cached.__init__')
        super().__init__(*args, **kwargs)
        self.__cache = weakref.WeakValueDictionary()

    def __call__(self, *args):
        # print('in Cached.__call__')
        if args in self.__cache:
            return self.__cache[args]
        else:
            obj = super().__call__(*args)
            self.__cache[args] = obj
            return obj


def make_sig(*names):
    """Helper function for the `ClassSignatureMeta` class."""
    params = [Parameter(name, Parameter.POSITIONAL_OR_KEYWORD)
              for name in names]
    return Signature(params)


class ClassSignatureMeta(type):
    def __new__(cls, clsname, bases, clsdict):
        clsdict['__signature__'] = make_sig(*clsdict.get('_fields', []))
        return super().__new__(cls, clsname, bases, clsdict)


class ClassSignature(metaclass=ClassSignatureMeta):
    _fields = []

    def __init__(self, *args, **kwargs):
        bound_values = self.__signature__.bind(*args, **kwargs)
        for name, value in bound_values.arguments.items():
            setattr(self, name, value)


class NoInstances(type):
    """Not callable class type."""
    def __call__(self, *args, **kwargs):
        raise RuntimeError("Can't instantiate directly.")


class Singleton(type):
    """Singleton class type."""
    def __init__(self, *args, **kwargs):
        # print('in Singleton.__init__')
        self.__instance = None
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        # print('in Singleton.__call__')
        if self.__instance is None:
            self.__instance = super().__call__(*args, **kwargs)
            return self.__instance
        else:
            return self.__instance


class method_func:
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
        sig = self.__name__ + get_object_signature(meth)
        if meth is not None:
            doc = """    {}\n{}""".format(sig,
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


class with_doc:
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
