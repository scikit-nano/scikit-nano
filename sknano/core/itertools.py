# -*- coding: utf-8 -*-
"""
===============================================================================
Custom iterator functions (:mod:`sknano.core.itertools`)
===============================================================================

.. currentmodule:: sknano.core.itertools

The functions defined in this module are from the
|itertools-recipes|_ for Python 3 or from the
Python Cookbook, Third Edition [PythonCookbook]_.

.. |itertools-recipes| replace:: Itertools Recipes
.. _itertools-recipes:
   https://docs.python.org/3/library/itertools.html#itertools-recipes
.. [PythonCookbook] Python Cookbook, Third Edition
   by David Beazley and Brian K. Jones.

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from collections import Iterable
from operator import itemgetter

import operator
import collections
import random
from itertools import chain, combinations, count, cycle, filterfalse, \
    groupby, islice, repeat, starmap, tee, zip_longest

__all__ = ['cyclic_pairs', 'take', 'tabulate', 'tail', 'consume', 'nth',
           'all_equal', 'quantify', 'padnone', 'ncycles', 'dotproduct',
           'flatten', 'repeatfunc', 'pairwise', 'grouper', 'roundrobin',
           'partition', 'powerset', 'unique_elements', 'unique_everseen',
           'unique_justseen', 'iter_except', 'first_true', 'random_product',
           'random_permutation', 'random_combination',
           'random_combination_with_replacement']


def cyclic_pairs(iterable):
    """Return iterator over all subsequence pairs in `iterable`, up to
    last and first element pair.

    Returns
    -------
    :class:`~python:collections.abc.Iterator`

    Examples
    --------
    >>> list(cyclic_pairs('ABC'))
    [('A', 'B'), ('B', 'C'), ('C', 'A')]

    """
    a, b = tee(iterable)
    return zip(a, chain(b, [next(b)]))


def take(n, iterable):
    """Return first n items of the iterable as a list.

    Parameters
    ----------
    n : :class:`~python:int`
    iterable : :class:`collections.abc.Iterable`

    Returns
    -------
    :class:`~python:list`

    Examples
    --------
    >>> from sknano.core import take, tabulate
    >>> t = tabulate(lambda i: i)
    >>> take(5, t)
    [0, 1, 2, 3, 4]
    >>> take(5, t)
    [5, 6, 7, 8, 9]

    """
    return list(islice(iterable, n))


def tabulate(function, start=0):
    """Return an iterator mapping `function` over linear input.

    Returns function(0), function(1), ...

    Parameters
    ----------
    function : callable
    start : :class:`~python:int`
        The start argument will be increased by 1 each time the iterator
        is called and fed into `function`.

    Returns
    -------
    :class:`~python:collections.abc.Iterator`

    Examples
    --------
    >>> from sknano.core import tabulate, ordinal_form, take
    >>> t = tabulate(ordinal_form, 10)
    >>> take(5, t)
    ['10th', '11th', '12th', '13th', '14th']
    >>> take(5, t)
    ['15th', '16th', '17th', '18th', '19th']
    >>> t = tabulate(ordinal_form)
    >>> take(5, t)
    ['0th', '1st', '2nd', '3rd', '4th']
    >>> take(5, t)
    ['5th', '6th', '7th', '8th', '9th']

    """
    return map(function, count(start))


def tail(n, iterable):
    """Return an iterator over the last `n` items.

    Parameters
    ----------
    n : :class:`~python:int`
    iterable : :class:`~python:collections.abc.Iterable`

    Returns
    -------
    :class:`~python:collections.abc.Iterator`.

    """
    return iter(collections.deque(iterable, maxlen=n))


def consume(iterator, n=None):
    """Advance the iterator n-steps ahead. If n is `None`, consume entirely.

    Parameters
    ----------
    iterator : :class:`~python:collections.abc.Iterator`
    n : {:class:`~python:None`, :class:`~python:int`}, optional

    """
    if n is None:
        # feed the entire iterator into a zero-length deque
        collections.deque(iterator, maxlen=0)
    else:
        # advance to the empty slice starting at position n
        next(islice(iterator, n, n), None)


def nth(iterable, n, default=None):
    """Returns the nth item or a default value.

    Parameters
    ----------
    iterable : :class:`~python:collections.abc.Iterable`
    n : :class:`~python:int`
    default : {:class:`~python:None`, :class:`~python:object`}, optional

    Examples
    --------
    >>> from sknano.core import nth
    >>> nth(range(10), 3)
    3
    >>> nth(range(10), 20, default='donkey')
    'donkey'

    """
    return next(islice(iterable, n, None), default)


def all_equal(iterable):
    """Returns True if all the elements are equal to each other.

    Parameters
    ----------
    iterable : :class:`~python:collections.abc.Iterable`

    Returns
    -------
    :class:`~python:bool`

    """
    g = groupby(iterable)
    return next(g, True) and not next(g, False)


def quantify(iterable, pred=bool):
    """Count how many times the predicate is true.

    Examples
    --------
    >>> from sknano.core import quantify
    >>> quantify([True, False, True])
    2

    """
    return sum(map(pred, iterable))


def padnone(iterable):
    """Returns the sequence elements and then returns None indefinitely.

    Examples
    --------
    >>> from sknano.core import padnone
    >>> take(5, padnone(range(3)))
    [0, 1, 2, None, None]

    """
    return chain(iterable, repeat(None))


def ncycles(iterable, n):
    """Returns the sequence elements n times"""
    return chain.from_iterable(repeat(tuple(iterable), n))


def dotproduct(vec1, vec2):
    """Returns the dot product of two \
        :class:`~python:collections.abc.Iterable`\ s.

    Parameters
    ----------
    vec1, vec2 : :class:`~python:collections.abc.Iterable`

    """
    return sum(map(operator.mul, vec1, vec2))


def flatten(items, ignore_types=(str, bytes)):
    """Yields single sequence of values with no nesting :ref:`[PythonCookbook]`

    Parameters
    ----------
    items : sequence
        nested sequence that you want to flatten into a single sequence with
        no nesting
    ignore_types : :class:`~python:tuple`, optional
        :class:`~python:tuple` of :class:`~python:collections.abc.Iterable`
        :class:`~python:type`\ s to that should not be interpreted as
        :class:`~python:collections.abc.Iterable`\ s to be flattened.
        Default is (:class:`~python:str`, :class:`~python:bytes`) to
        prevent strings and bytes from being interpreted as iterables
        and expanded as individual characters.

    Yields
    ------
    sequence
        single sequence with no nesting

    Examples
    --------
    >>> items = [1, 2, [3, 4, [5, 6], 7], 8]
    >>> from sknano.core import flatten
    >>> items = list(flatten(items))
    >>> print(items)
    [1, 2, 3, 4, 5, 6, 7, 8]

    """
    for x in items:
        if isinstance(x, Iterable) and not isinstance(x, ignore_types):
            yield from flatten(x)
        else:
            yield x


def repeatfunc(func, times=None, *args):
    """Repeat calls to func with specified arguments.

    """
    if times is None:
        return starmap(func, repeat(args))
    return starmap(func, repeat(args, times))


def pairwise(iterable):
    """Returns an iterator of paired items, overlapping, from the original.

    Examples
    --------
    >>> from sknano.core import pairwise, tabulate, take
    >>> t = tabulate(lambda i: i)
    >>> take(5, pairwise(t))
    [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]
    >>> take(5, pairwise(t))
    [(6, 7), (7, 8), (8, 9), (9, 10), (10, 11)]

    """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks.

    Examples
    --------
    >>> from sknano.core import grouper
    >>> list(grouper('ABCDEFG', 3, '?'))
    [('A', 'B', 'C'), ('D', 'E', 'F'), ('G', '?', '?')]

    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def roundrobin(*iterables):
    """Yields an item from each iterable, alternating between them.

    Recipe credited to George Sakkis

    Examples
    --------
    >>> from sknano.core import roundrobin
    >>> list(roundrobin('ABC', 'D', 'EF'))
    ['A', 'D', 'E', 'B', 'F', 'C']

    """
    pending = len(iterables)
    nexts = cycle(iter(it).__next__ for it in iterables)
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1
            nexts = cycle(islice(nexts, pending))


def partition(pred, iterable):
    """Use a predicate to partition entries into false entries and true \
    entries.

    Parameters
    ----------
    pred : callable
        callable function taking one argument and returning
        :class:`~python:True` or :class:`~python:False`.
    iterable : :class:`~python:collections.abc.Iterable`

    Returns
    -------
    :class:`~python:tuple` of :class:`~pyton:collections.abc.Iterables`

    Examples
    --------
    >>> from sknano.core import partition
    >>> is_odd = lambda x: x % 2 != 0
    >>> p = partition(is_odd, range(10))
    >>> [print(i) for i in p]
    [[0, 2, 4, 6, 8], [1, 3, 5, 7, 9]]

    """
    t1, t2 = tee(iterable)
    return filterfalse(pred, t1), filter(pred, t2)


def powerset(iterable):
    """Yields all possible subsets of the iterable.

    Examples
    --------
    >>> from sknano.core import powerset
    >>> list(powerset([1, 2, 3]))
    [(), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]

    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def unique_everseen(iterable, key=None):
    """Yields unique elements in `iterable`, preserving order, remembering \
        all elements ever seen.

    Parameters
    ----------
    iterable : :class:`~python:collections.abc.Iterable`
    key : {None, callable}, optional

    Yields
    ------
    element in iterable

    Examples
    --------
    >>> list(unique_everseen('AAAABBBCCDAABBB'))
    ['A', 'B', 'C', 'D']
    >>> list(unique_everseen('ABBCcAD'))
    ['A', 'B', 'C', 'c', 'D']
    >>> list(unique_everseen('ABBCcAD', str.lower))
    ['A', 'B', 'C', 'D']

    """
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in filterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element


def unique_elements(iterable, key=None):
    """Alias for :func:`~sknano.core.unique_everseen`."""
    return unique_everseen(iterable, key)


def unique_justseen(iterable, key=None):
    """Yields unique elements in `iterable`, preserving order, remembering \
        only the element just seen.

    Parameters
    ----------
    iterable : :class:`~python:collections.abc.Iterable`
    key : {None, callable}, optional

    Yields
    ------
    element in iterable

    Examples
    --------
    >>> list(unique_justseen('AAAABBBCCDAABBB'))
    ['A', 'B', 'C', 'D', 'A', 'B']
    >>> list(unique_justseen('ABBCcAD', str.lower))
    ['A', 'B', 'C', 'A', 'D']

    """
    return map(next, map(itemgetter(1), groupby(iterable, key)))


def iter_except(func, exception, first=None):
    """Call a function repeatedly until an exception is raised.

    Converts a call-until-exception interface to an iterator interface.
    Like builtins.iter(func, sentinel) but uses an exception instead
    of a sentinel to end the loop.

    Parameters
    ----------
    func : `callable`
    exception : :class:`~python:Exception`
    first : `callable`
        for database APIs needing an initial cast to db.first()

    Examples
    --------

    >>> from sknano.core import iter_except
    >>> list(iter_except(list(range(3)).pop, IndexError))
    [2, 1, 0]
    >>>

    priority queue iterator:

    >>> # iter_except(functools.partial(heappop, h), IndexError)

    non-blocking dict iterator:

    >>> # iter_except(d.popitem, KeyError)

    non-blocking deque iterator:

    >>> # iter_except(d.popleft, IndexError)

    loop over a producer Queue:

    >>> # iter_except(q.get_nowait, Queue.Empty)

    non-blocking set iterator:

    >>> # iter_except(s.pop, KeyError)

    """
    try:
        # For database APIs needing an initial cast to db.first()
        if first is not None:
            yield first()
        while 1:
            yield func()
    except exception:
        pass


def first_true(iterable, default=False, pred=None):
    """Returns the first true value in the iterable.

    If no true value is found, returns *default*

    If *pred* is not None, returns the first item
    for which pred(item) is true.

    Parameters
    ----------
    default : `bool`
    pred : {None, `callable`}, optional

    Examples
    --------
    >>> from sknano.core import first_true
    >>> first_true([(), None, '', 0, 1])
    1
    >>> first_true([(), None, '', 0], 'donkey')
    'donkey'
    >>> is_odd = lambda x: x % 2 != 0
    >>> first_true([2, 4, 6, 8, 11], default='donkey', pred=is_odd)
    11
    >>> first_true([2, 4, 6, 8], default='donkey', pred=is_odd)
    'donkey'
    """
    return next(filter(pred, iterable), default)


def random_product(*args, repeat=1):
    """Random selection from :func:`~python:itertools.product`."""
    pools = [tuple(pool) for pool in args] * repeat
    return tuple(random.choice(pool) for pool in pools)


def random_permutation(iterable, r=None):
    """Random selection from :func:`~python:itertools.permutations`."""
    pool = tuple(iterable)
    r = len(pool) if r is None else r
    return tuple(random.sample(pool, r))


def random_combination(iterable, r):
    """Random selection from :func:`~python:itertools.combinations`."""
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return tuple(pool[i] for i in indices)


def random_combination_with_replacement(iterable, r):
    """Random selection from \
    :func:`~python:itertools.combinations_with_replacement`.

    """
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.randrange(n) for i in range(r))
    return tuple(pool[i] for i in indices)
