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

__all__ = ['plural_word_check', 'pluralize', 'ordinal_form']


def plural_word_check(word, count):
    """Alias for :func:`pluralize`"""
    return pluralize(word, count)


def pluralize(word, count):
    """Make a word plural by adding an *s* if `count` != 1.

    Parameters
    ----------
    word : :class:`~python:str`
        the word
    count : :class:`~python:int`
        the word count

    Returns
    -------
    :class:`~python:str`

    """
    return word if count == 1 else word + 's'


def ordinal_form(n):
    """Convert number to ordinal form in English.

    Parameters
    ----------
    n : :class:`~python:int`

    Returns
    -------
    :class:`~python:str`

    Examples
    --------
    >>> from pprint import pprint
    >>> from sknano.core import ordinal_form
    >>> pprint([ordinal_form(i) for i in range(200)], width=70, compact=True)
    ['0th', '1st', '2nd', '3rd', '4th', '5th', '6th', '7th', '8th', '9th',
     '10th', '11th', '12th', '13th', '14th', '15th', '16th', '17th',
     '18th', '19th', '20th', '21st', '22nd', '23rd', '24th', '25th',
     '26th', '27th', '28th', '29th', '30th', '31st', '32nd', '33rd',
     '34th', '35th', '36th', '37th', '38th', '39th', '40th', '41st',
     '42nd', '43rd', '44th', '45th', '46th', '47th', '48th', '49th',
     '50th', '51st', '52nd', '53rd', '54th', '55th', '56th', '57th',
     '58th', '59th', '60th', '61st', '62nd', '63rd', '64th', '65th',
     '66th', '67th', '68th', '69th', '70th', '71st', '72nd', '73rd',
     '74th', '75th', '76th', '77th', '78th', '79th', '80th', '81st',
     '82nd', '83rd', '84th', '85th', '86th', '87th', '88th', '89th',
     '90th', '91st', '92nd', '93rd', '94th', '95th', '96th', '97th',
     '98th', '99th', '100th', '101st', '102nd', '103rd', '104th', '105th',
     '106th', '107th', '108th', '109th', '110th', '111st', '112nd',
     '113rd', '114th', '115th', '116th', '117th', '118th', '119th',
     '120th', '121st', '122nd', '123rd', '124th', '125th', '126th',
     '127th', '128th', '129th', '130th', '131st', '132nd', '133rd',
     '134th', '135th', '136th', '137th', '138th', '139th', '140th',
     '141st', '142nd', '143rd', '144th', '145th', '146th', '147th',
     '148th', '149th', '150th', '151st', '152nd', '153rd', '154th',
     '155th', '156th', '157th', '158th', '159th', '160th', '161st',
     '162nd', '163rd', '164th', '165th', '166th', '167th', '168th',
     '169th', '170th', '171st', '172nd', '173rd', '174th', '175th',
     '176th', '177th', '178th', '179th', '180th', '181st', '182nd',
     '183rd', '184th', '185th', '186th', '187th', '188th', '189th',
     '190th', '191st', '192nd', '193rd', '194th', '195th', '196th',
     '197th', '198th', '199th']

    """
    ordinal_suffix = {}
    ordinal_suffix.update(dict.fromkeys(range(20), 'th'))
    ordinal_suffix.update({1: 'st', 2: 'nd', 3: 'rd'})
    try:
        return ''.join((str(n), ordinal_suffix[n]))
    except KeyError:
        last_digit = int(str(n)[-1])
        return ''.join((str(n), ordinal_suffix[last_digit]))
