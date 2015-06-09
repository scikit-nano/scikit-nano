# -*- coding: utf-8 -*-
"""
===============================================================================
Custom iterator functions (:mod:`sknano.core._itertools`)
===============================================================================

.. currentmodule:: sknano.core._itertools

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import collections
import operator
import random
from itertools import chain, combinations, count, cycle, filterfalse, \
    islice, repeat, starmap, tee, zip_longest

__all__ = ['cyclic_pairs', 'take', 'tabulate', 'consume', 'nth', 'quantify',
           'padnone', 'ncycles', 'dotproduct', 'flatten', 'repeatfunc',
           'pairwise', 'grouper', 'roundrobin', 'partition', 'powerset',
           'unique_elements', 'iter_except', 'first_true',
           'random_product', 'random_permutation', 'random_combination',
           'random_combination_with_replacement']


def cyclic_pairs(iterable):
    """
    Generate a cyclic list of all subsequence pairs in `iterable`.

    Returns
    -------
    :class:`~python:list`

    Examples
    --------

    >>> cyclic_pairs('ABC')
    [('A', 'B'), ('B', 'C'), ('C', 'A')]

    """
    a, b = tee(iterable)
    return list(zip(a, chain(b, [next(b)])))


def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))


def tabulate(function, start=0):
    "Return function(0), function(1), ..."
    return map(function, count(start))


def consume(iterator, n):
    "Advance the iterator n-steps ahead. If n is none, consume entirely."
    # Use functions that consume iterators at C speed.
    if n is None:
        # feed the entire iterator into a zero-length deque
        collections.deque(iterator, maxlen=0)
    else:
        # advance to the empty slice starting at position n
        next(islice(iterator, n, n), None)


def nth(iterable, n, default=None):
    "Returns the nth item or a default value"
    return next(islice(iterable, n, None), default)


def quantify(iterable, pred=bool):
    "Count how many times the predicate is true"
    return sum(map(pred, iterable))


def padnone(iterable):
    """Returns the sequence elements and then returns None indefinitely.

    Useful for emulating the behavior of the built-in map() function.
    """
    return chain(iterable, repeat(None))


def ncycles(iterable, n):
    "Returns the sequence elements n times"
    return chain.from_iterable(repeat(tuple(iterable), n))


def dotproduct(vec1, vec2):
    return sum(map(operator.mul, vec1, vec2))


def flatten(listOfLists):
    "Flatten one level of nesting"
    return chain.from_iterable(listOfLists)


def repeatfunc(func, times=None, *args):
    """Repeat calls to func with specified arguments.

    Example:  repeatfunc(random.random)
    """
    if times is None:
        return starmap(func, repeat(args))
    return starmap(func, repeat(args, times))


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def roundrobin(*iterables):
    "roundrobin('ABC', 'D', 'EF') --> A D E B F C"
    # Recipe credited to George Sakkis
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
    "Use a predicate to partition entries into false entries and true entries"
    # partition(is_odd, range(10)) --> 0 2 4 6 8   and  1 3 5 7 9
    t1, t2 = tee(iterable)
    return filterfalse(pred, t1), filter(pred, t2)


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def unique_elements(iterable, key=None):
    """Return generator of unique elements in `iterable`, preserving order.

    Examples
    --------
    >>> list(unique_elements('AAAABBBCCDAABBB'))
    ['A', 'B', 'C', 'D']
    >>> list(unique_elements('ABBCcAD'))
    ['A', 'B', 'C', 'c', 'D']
    >>> list(unique_elements('ABBCcAD', str.lower))
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


def iter_except(func, exception, first=None):
    """Call a function repeatedly until an exception is raised.

    Converts a call-until-exception interface to an iterator interface.
    Like builtins.iter(func, sentinel) but uses an exception instead
    of a sentinel to end the loop.

    Parameters
    ----------
    func : `callable`
    exception : `Exception`
    first : `callable`
        for database APIs needing an initial cast to db.first()

    Examples
    --------
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
    >>> # first_true([a,b,c], x)
    >>> # a or b or c or x
    >>> # first_true([a,b], x, f)
    >>> # a if f(a) else b if f(b) else x

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
