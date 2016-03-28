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

from textwrap import shorten

from tabulate import tabulate

from pyparsing import Group, Forward, Optional, Regex, Suppress, Keyword, \
    Literal, Word, ZeroOrMore, alphas, alphanums, delimitedList, \
    oneOf, replaceWith, quotedString, removeQuotes

__all__ = ['TabulateMixin', 'pluralize', 'plural_word_check', 'ordinal_form',
           'map_operator', 'map_function', 'asint', 'asfloat', 'asbool',
           'astuple', 'aslist', 'asdict', 'asset', 'integer', 'real',
           'number', 'boolean', 'string', 'none', 'LPAR', 'RPAR',
           'LBRACK', 'RBRACK', 'LBRACE', 'RBRACE', 'COLON', 'SEMICOLON',
           'SPACE', 'COMMA', 'EQUAL', 'binary_operator', 'hashable_item',
           'unhashable_item', 'expr_item', 'tuple_expr', 'list_expr',
           'dict_expr', 'set_expr', 'kwarg', 'kwargs_expr', 'signature_args',
           'signature_kwargs', 'call_signature']


class TabulateMixin:
    """Mixin class for pretty tabulated output strings."""

    def __str__(self):
        title = self._table_title_str()
        table = self._tabulate()
        return '\n'.join((title, table))

    def _tabulate(self, values=None, headers=(), tablefmt='fancy_grid'):
        if values is None:
            return self._tabulate(*self._tabular_data())
        return tabulate(values, headers=headers, tablefmt=tablefmt)

    def _tabular_data(self):
        """Return :class:`~python:tuple` of tabular data for \
            pretty tabulated output."""
        print('in TabulateMixin._tabular_data()')
        header = self.__class__.__name__
        fmt = self._tabular_data_format_string
        return [fmt(self)], (header,)

    def _table_title_str(self, prepend='', obj=None, append='', sep=''):
        if obj is None:
            title = repr(self)
        elif obj is not None and not isinstance(obj, str):
            title = repr(obj)
        title = '{}'.format(shorten(title, width=78))
        return sep.join((prepend, title, append))

    def _tabular_data_format_string(self, value, begin=None, end=None):
        strrep = '{!r}'.format(value)
        try:
            return strrep[begin:end]
        except IndexError:
            return strrep


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

    Examples
    --------
    On occasion, it is desirable to describe a numerical value in terms of
    a noun which the number is quantifying. For example, given
    a function which accepts a numerical parameter `n` and returns
    a string describing the number of `n` *objects*, then this
    helper function may be of use. For example::

    >>> from sknano.core import pluralize
    >>> def apple_count(n):
    ...     return '{} {}'.format(n, pluralize('apple', n))
    ...
    >>> [apple_count(i) for i in range(3)]
    ['0 apples', '1 apple', '2 apples']

    """
    return word if count == 1 else word + 's'


def plural_word_check(word, count):
    """Alias for :func:`pluralize`"""
    return pluralize(word, count)


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

type_map = {
    'bool': bool,
    'int': int,
    'float': float,
    'list': list,
    'str': str,
    'tuple': tuple,
}


def map_type_str(type_str):
    return type_map[type_str]


def cast_type(parameters):
    pass


def asint(s, l, t):
    """:mod:`~pyparsing:pyparsing` parse action to convert token to \
        :class:`~python:int`."""
    return int(t[0])


def asfloat(s, l, t):
    """:mod:`~pyparsing:pyparsing` parse action to convert token to \
        :class:`~python:float`."""
    try:
        val = int(t[0])
    except ValueError:
        val = float(t[0])
    return val


def asbool(s, l, t):
    """:mod:`~pyparsing:pyparsing` parse action to convert token to \
        :class:`~python:bool`."""
    return t[0] == 'True'


def astuple(s, l, t):
    """:mod:`~pyparsing:pyparsing` parse action to convert token to \
        :class:`~python:tuple`."""
    return tuple(t.asList())


def aslist(s, l, t):
    """:mod:`~pyparsing:pyparsing` parse action to convert token to \
        :class:`~python:list`."""
    return [t.asList()]


def asdict(s, l, t):
    """:mod:`~pyparsing:pyparsing` parse action to convert token to \
        :class:`~python:dict`."""
    return dict(t.asList())


def asset(s, l, t):
    """:mod:`~pyparsing:pyparsing` parse action to convert token to \
        :class:`~python:set`."""
    return set(t.asList())


def map_function(s, l, t):
    """:mod:`~pyparsing:pyparsing` parse action to convert function string to \
        function in :func:`sknano.core.math.function_map`."""
    from sknano.core.math import function_map
    return function_map[t[0]]


def map_operator(s, l, t):
    """:mod:`~pyparsing:pyparsing` parse action to convert operator string to \
        operator in :func:`sknano.core.math.operator_map`."""
    from sknano.core.math import operator_map
    return operator_map[t[0]]


integer = Regex(r'[+-]?\d+').setParseAction(asint)
real = Regex(r'[+-]?\d+(\.\d*)?([eE][+-]?\d+)?').setParseAction(asfloat)
number = real | integer

boolean = oneOf("True False").setParseAction(asbool)
none = Literal("None").setParseAction(replaceWith(None))
string = quotedString.setParseAction(removeQuotes)

# point = Literal('.')
# e = CaselessLiteral('E')
# plus_or_minus_sign = Word('+-', exact=1)
# number = Word(nums)
# positive_integer = number.setParseAction(asint)
# real_number = \
#     Combine(Optional(plus_or_minus_sign) +
#             (number + point + Optional(number) | (point + number)) +
#             Optional(e) + Optional(plus_or_minus_sign) + number) \
#     .setParseAction(asfloat)

LPAR, RPAR = map(Suppress, '()')
LBRACK, RBRACK = map(Suppress, '[]')
LBRACE, RBRACE = map(Suppress, '{}')
COLON = Suppress(':')
SEMICOLON = Suppress(';')
SPACE = Suppress(' ')
COMMA = Suppress(',')
EQUAL = Suppress('=')

binary_operator = oneOf('< <= == > >= != LT LE EQ GT GE NE', caseless=True)
binary_operator.setParseAction(map_operator)

tuple_expr = Forward()
list_expr = Forward()
dict_expr = Forward()
set_expr = Forward()

hashable_item = real | integer | string | boolean | none | tuple_expr
unhashable_item = list_expr | dict_expr | set_expr
expr_item = hashable_item | unhashable_item

tuple_expr << (LPAR +
               Optional(delimitedList(expr_item) + Optional(COMMA)) +
               RPAR).setParseAction(astuple)

list_expr << (LBRACK +
              Optional(delimitedList(expr_item) + Optional(COMMA)) +
              RBRACK).setParseAction(aslist)

kwarg = Group(Word(alphas, alphanums + '_') + EQUAL + expr_item)
kwargs_expr = Forward()
kwargs_expr << delimitedList(kwarg).setParseAction(asdict)

dict_entry = Group(hashable_item + COLON + expr_item) | kwarg
dict_entries = Optional(delimitedList(dict_entry) + Optional(COMMA))
dict_expr << (LBRACE + dict_entries + RBRACE |
              Suppress(Keyword('dict')) + LPAR +
              dict_entries + RPAR).setParseAction(asdict)


signature_args = Forward()
signature_args << delimitedList(expr_item).setParseAction(aslist)
signature_kwargs = kwargs_expr.copy()

call_signature = Group(Optional(signature_args + COMMA) +
                       ZeroOrMore(signature_kwargs))
