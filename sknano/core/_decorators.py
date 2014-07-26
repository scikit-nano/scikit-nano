# -*- coding: utf-8 -*-
"""
======================================================================
Decorator classes (:mod:`sknano.core._decorators`)
======================================================================

.. currentmodule:: sknano.core._decorators

Contents
========

"""
from __future__ import absolute_import, division, print_function

__all__ = ['deprecated', 'with_doc']


import warnings
import functools


def deprecated(func):
    """Deprecated decorator.

    A decorator which can be used to mark functions as deprecated.
    replacement is a callable that will be called with the same args
    as the decorated function.

    Source: http://code.activestate.com/recipes/577819-deprecated-decorator/
    Original Author: Giampaolo Rodola' <g.rodola [AT] gmail [DOT] com>
    License: MIT

    Returns
    -------
    wrapper : decorator function

    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        warnings.warn_explicit(
            "Call to deprecated function {}".format(func.__name__),
            category=DeprecationWarning,
            filename=func.func_code.co_filename,
            lineno=func.func_code.co_firstlineno + 1
        )
        return func(*args, **kwargs)
    return wrapper


def deprecated_class(replacement=None):
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
    def decorator(func):
        msg = "{!s} is deprecated".format(func.__name__)
        if replacement is not None:
            msg += "; use {!s} instead".format(replacement)
        if func.__doc__ is None:
            func.__doc__ = msg

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            warnings.warn(msg, category=DeprecationWarning, stacklevel=2)
            return func(*args, **kwargs)
        return wrapper
    return decorator


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
