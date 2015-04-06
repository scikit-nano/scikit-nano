# -*- coding: utf-8 -*-
"""
=========================================================================
Helper functions for string manipulation (:mod:`sknano.core._strings`)
=========================================================================

.. currentmodule:: sknano.core._strings

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

__all__ = ['plural_word_check', 'pluralize']


def plural_word_check(word, count):
    """Alias for :func:`pluralize`"""
    return pluralize(word, count)


def pluralize(word, count):
    """Make a word plural by adding an *s* if `count` != 1.

    Parameters
    ----------
    word : str
        the word
    count : int
        the word count

    Returns
    -------
    str

    """
    return word if count == 1 else word + 's'
