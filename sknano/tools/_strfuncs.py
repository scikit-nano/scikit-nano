# -*- coding: utf-8 -*-
"""
=========================================================================
Helper functions for string manipulation (:mod:`sknano.tools._strfuncs`)
=========================================================================

.. currentmodule:: sknano.tools._strfuncs

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

__all__ = ['plural_word_check']


def plural_word_check(word, count):
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
