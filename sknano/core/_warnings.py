# -*- coding: utf-8 -*-
"""
======================================================================
Decorator classes (:mod:`sknano.core._warnings`)
======================================================================

.. currentmodule:: sknano.core._warnings

Contents
========

"""
from __future__ import absolute_import, division, print_function

__all__ = ['package_removed_warning']


import warnings

def package_removed_warning(oldpkg, newpkg):
    return warnings.warn(
        "The {!r} package was removed or renamed in ".format(oldpkg) +
        "version 0.3.\n"
        "Most modules previously import from {!r},\n".format(oldpkg) +
        "can now be found in the {!r} package.".format(newpkg) +
        "Please check the documentation and review the module specific\n"
        "ImportWarnings and update your code accordingly.\n"
        "The fallback code and package ImportWarnings will be completely\n"
        "removed in a future version.\n",
        ImportWarning
    )
