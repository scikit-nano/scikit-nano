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

__all__ = ['removed_package_warning']


import warnings

warnings.resetwarnings()


def removed_package_warning(oldpkg, newpkg=None):
    msg = '\n\n#TODO: CREATE CUSTOM TEMPLATE STRING FOR THIS WARNING MESSAGE.\n'
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
