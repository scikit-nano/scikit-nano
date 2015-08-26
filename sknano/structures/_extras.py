# -*- coding: utf-8 -*-
"""
===============================================================================
Extra helper functions (:mod:`sknano.structures._extras`)
===============================================================================

.. currentmodule:: sknano.structures._extras

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import importlib
import numbers
import re

import numpy as np

from sknano.core.math import comparison_symbol_operator_mappings

__all__ = ['cmp_Ch', 'filter_Ch', 'filter_Ch_list', 'generate_Ch_list',
           'generate_Ch_property_grid', 'get_Ch_indices',
           'get_chiral_type', 'get_Ch_type',
           'map_Ch', 'chiral_type_name_mappings', 'CHIRAL_TYPES',
           'filter_key_type_mappings', 'attr_units', 'attr_symbols',
           'attr_strfmt', 'type_check_chiral_indices',
           'get_chiral_indices', 'get_chiral_indices_from_str']

chiral_type_name_mappings = \
    {'achiral': 'aCh', 'armchair': 'AC', 'zigzag': 'ZZ', 'chiral': 'Ch'}
CHIRAL_TYPES = list(chiral_type_name_mappings.keys()) + \
    list(chiral_type_name_mappings.values())

filter_key_type_mappings = {}
filter_key_type_mappings['Ch_type'] = str
for k in ('even', 'odd'):
    filter_key_type_mappings[k + '_only'] = bool

for k in ('min_index', 'max_index',
          'min_n', 'max_n',
          'min_m', 'max_m',
          'n', 'm'):
    filter_key_type_mappings[k] = int

attr_units = {}
attr_units['dt'] = \
    attr_units['rt'] = \
    attr_units['Ch'] = \
    attr_units['T'] = \
    attr_units['bond'] = ' \u212B'
attr_units['chiral_angle'] = '\u00b0'

attr_symbols = {}
attr_symbols['t1'] = 't\u2081'
attr_symbols['t2'] = 't\u2082'
attr_symbols['chiral_angle'] = '\u03b8c'

attr_strfmt = {}
attr_strfmt['Ch'] = \
    attr_strfmt['T'] = \
    attr_strfmt['dt'] = \
    attr_strfmt['rt'] = \
    attr_strfmt['chiral_angle'] = '{:.2f}'
attr_strfmt['bond'] = '{:.3f}'


def cmp_Ch(Ch1, Ch2):
    """Custom comparator function for sorting chirality lists.

    Parameters
    ----------
    Ch1, Ch2 : {str, tuple}
        2-tuple or 2-tuple string of integers

    Returns
    -------
    int

    """

    if isinstance(Ch1, tuple) and isinstance(Ch2, tuple):
        n1, m1 = Ch1
        Ch1_type = get_Ch_type(Ch1)
        n2, m2 = Ch2
        Ch2_type = get_Ch_type(Ch2)
    elif isinstance(Ch1, str) and \
            isinstance(Ch2, str):
        n1, m1 = get_chiral_indices_from_str(Ch1)
        Ch1_type = get_Ch_type(Ch1)
        n2, m2 = get_chiral_indices_from_str(Ch2)
        Ch2_type = get_Ch_type(Ch2)

    if Ch1_type == 'armchair' and Ch2_type == 'zigzag':
        return 1
    elif Ch1_type == 'zigzag' and Ch2_type == 'armchair':
        return -1
    else:
        if n1 > n2:
            return 1
        elif n1 < n2:
            return -1
        else:
            if m1 > m2:
                return 1
            elif m1 < m2:
                return -1
            else:
                return 0


def filter_Ch(Ch, even_only=False, odd_only=False, Ch_type=None,
              min_index=None, max_index=None, min_n=None, max_n=None,
              min_m=None, max_m=None, **kwargs):
    """Filter for testing if chirality satisfies given constraint parameters.

    Parameters
    ----------
    Ch : {str, tuple}
    even_only, odd_only : bool, optional
    Ch_type : {None, 'armchair', 'zigzag', 'achiral', 'chiral'}, optional
    min_index, max_index : int, optional
    min_n, max_n : int, optional
    min_m, max_m : int, optional

    Returns
    -------
    bool

    """

    if isinstance(Ch, tuple):
        this_Ch_type = get_Ch_type(Ch)
        n, m = Ch
    else:
        n, m = get_chiral_indices_from_str(Ch)
        this_Ch_type = get_Ch_type(Ch)

    if even_only:
        if n % 2 == 0 and m % 2 == 0:
            return True
        else:
            return False

    elif odd_only:
        if this_Ch_type != 'zigzag':
            if n % 2 != 0 and m % 2 != 0:
                return True
            else:
                return False
        else:
            if max(n, m) % 2 != 0:
                return True
            else:
                return False

    elif Ch_type is not None:
        if this_Ch_type in ('armchair', 'zigzag'):
            if Ch_type in ('achiral', 'aCh') or \
                    (Ch_type in ('armchair', 'AC') and
                     this_Ch_type == 'armchair') or \
                    (Ch_type in ('zigzag', 'ZZ') and
                     this_Ch_type == 'zigzag'):
                return True
            else:
                return False
        elif Ch_type in ('chiral', 'Ch') and this_Ch_type == 'chiral':
            return True
        else:
            return False

    elif min_index is not None and isinstance(min_index, int):
        if this_Ch_type != 'zigzag':
            if n >= min_index and m >= min_index:
                return True
            else:
                return False
        else:
            if max(n, m) >= min_index:
                return True
            else:
                return False

    elif max_index is not None and isinstance(max_index, int):
        if n <= max_index and m <= max_index:
            return True
        else:
            return False

    else:
        return True


def filter_Ch_list(Ch_list, property_filters=None, **kwargs):
    """Filter list of chiralities.

    Parameters
    ----------
    Ch_list : sequence
        list of chiralities
    property_filters : sequence, optional
    kwargs : dict, optional
        dictionary of conditions to test each chirality against.

    Returns
    -------
    :class:`python:list`
        filtered `Ch_list`

    """
    # Ch_list = np.asarray(Ch_list)
    if property_filters is not None:
        from ._swnt import SWNT
        filtered_list = Ch_list[:]
        try:
            for filter_index, (prop, cmp_symbol, value) in \
                    enumerate(property_filters, start=1):
                cmp_op = comparison_symbol_operator_mappings[cmp_symbol]
                tmp_list = []
                for Ch in filtered_list:
                    n, m = Ch
                    nanotube = SWNT(n=n, m=m)
                    try:
                        if cmp_op(getattr(nanotube, prop), value):
                            tmp_list.append(Ch)
                    except AttributeError:
                        break
                filtered_list = tmp_list[:]
        except ValueError as e:
            print(e)
        finally:
            Ch_list = filtered_list[:]
    return [Ch for Ch in Ch_list if filter_Ch(Ch, **kwargs)]


def generate_Ch_list(ns=None, ni=None, nf=None, dn=None,
                     ms=None, mi=None, mf=None, dm=None,
                     imax=None, chiral_types=None, handedness=None,
                     echo_zsh_str=False):
    """Generate a list of :math:`(n, m)` chiralities.

    Parameters
    ----------
    ns, ms : sequence
        list of :math:`n` and :math:`m` chiral indices.
    ni, nf, dn : int, optional
        :math:`(n_i, n_f, \\Delta n)` denote the `start`, `stop`, and
        `step` parameters passed to the numpy function `np.arange`,
        to generate an array of evenly spaced :math:`n` chiral indices.
        `ni` only required if `ns` sequence kwarg is `None`.
    mi, mf, dm : int, optional
        :math:`(m_i, m_f, \\Delta m)` denote the `start`, `stop`, and
        `step` parameters passed to the numpy function `np.arange`,
        to generate an array of evenly spaced :math:`m` chiral indices.
        `mi` only required if `ms` sequence kwarg is `None`.
    imax : int, optional
        maximum chiral index :math:`n = m = i_{max}`.
    chiral_types : {'all', 'achiral', 'chiral', 'armchair', 'zigzag'}, optional
    handedness : {'all', 'left', 'right'}, optional
    echo_zsh_str : bool, optional

    Returns
    -------
    Ch_list : sequence
        list of chiralities

    """
    Ch_list = []
    if imax is not None:
        nf = mf = imax

    try:
        if (not ns or (isinstance(ns, list) and ns[0] is None)) \
                and not ni and not nf and \
                (not ms or (isinstance(ns, list) and ms[0] is None)) \
                and not mi and not mf:
            raise TypeError('ns or [ni,] nf[, dn,] and/or '
                            'ms or [mi,] mf[, dm,] must be specified')
        if not ns and isinstance(ms, list) and ms[0] is not None and \
                not ni and not nf and not mi and not mf:
            ns = ms[:]
            ms = None
        elif not ns and (not ms or isinstance(ms, list) and ms[0] is None) \
                and not ni and not nf and \
                (isinstance(mi, int) or isinstance(mf, int)):
            ni = mi
            nf = mf
            mi = mf = None
            if not dn and isinstance(dm, int):
                dn = dm
                dm = None
        if not ns and (isinstance(ni, int) or isinstance(nf, int)):
            try:
                ns = np.arange(ni, nf + 1, dn, dtype=int).tolist()
            except TypeError:
                try:
                    ni, nf, dn = int(float(ni)), int(float(nf)), int(float(dn))
                    ns = np.arange(ni, nf + 1, dn, dtype=int).tolist()
                except (TypeError, ValueError):
                    try:
                        ni, nf = int(float(ni)), int(float(nf))
                        ns = np.arange(ni, nf + 1, dtype=int).tolist()
                    except (TypeError, ValueError):
                        try:
                            nf = int(float(nf))
                            ns = np.arange(nf + 1, dtype=int).tolist()
                        except (TypeError, ValueError):
                            ns = [int(float(ni))]
        if (not ms or isinstance(ms, list) and ms[0] is None) and \
                not mi and not mf:
            for n in ns:
                m = 0
                while m <= n:
                    if (m == 0) or (m == n):
                        Ch_list.append((n, m))
                    else:
                        if handedness == 'left':
                            Ch_list.append((m, n))
                        elif handedness == 'right':
                            Ch_list.append((n, m))
                        else:
                            Ch_list.append((n, m))
                            Ch_list.append((m, n))
                    m += 1
        elif (not ms or isinstance(ms, list) and ms[0] is None) and \
                (isinstance(mi, int) or isinstance(mf, int)):
            try:
                ms = np.arange(mi, mf + 1, dm, dtype=int).tolist()
            except TypeError:
                try:
                    mi, mf, dm = int(float(mi)), int(float(mf)), int(float(dm))
                    ms = np.arange(mi, mf + 1, dm, dtype=int).tolist()
                except (TypeError, ValueError):
                    try:
                        mi, mf = int(float(mi)), int(float(mf))
                        ms = np.arange(mi, mf + 1, dtype=int).tolist()
                    except (TypeError, ValueError):
                        try:
                            mf = int(float(mf))
                            ms = np.arange(mf + 1, dtype=int).tolist()
                        except (TypeError, ValueError):
                            ms = [int(float(mi))]
    except (TypeError, ValueError):
        Ch_list = None
    else:
        if len(Ch_list) == 0 and isinstance(ms, list) and len(ms) != 0:
            for n in ns:
                for m in ms:
                    if (n == 0) and (m == 0):
                        break
                    elif (n == 0) or (m == 0):
                        nm = (n, m)
                        mn = (m, n)
                        if handedness == 'left':
                            if m != 0 and nm not in Ch_list:
                                Ch_list.append(nm)
                        elif handedness == 'right':
                            if n != 0 and nm not in Ch_list:
                                Ch_list.append(nm)
                        else:
                            if nm not in Ch_list:
                                Ch_list.append(nm)
                            if mn not in Ch_list:
                                Ch_list.append(mn)
                    elif (n == m) and (n, m) not in Ch_list:
                        Ch_list.append((n, m))
                    else:
                        nm = (n, m)
                        mn = (m, n)
                        if handedness == 'left':
                            if n < m and nm not in Ch_list:
                                Ch_list.append(nm)
                        elif handedness == 'right':
                            if n > m and nm not in Ch_list:
                                Ch_list.append(nm)
                        else:
                            if nm not in Ch_list:
                                Ch_list.append(nm)
                            if mn not in Ch_list:
                                Ch_list.append(mn)
    finally:
        if echo_zsh_str and Ch_list is not None:
            print(' '.join([repr(str(Ch)) for Ch in Ch_list]))
        if chiral_types is not None:
            if not isinstance(chiral_types, list):
                chiral_types = [chiral_types]
            chiral_types = chiral_types[:]

            if any([chiral_type in CHIRAL_TYPES for chiral_type in
                    chiral_types]):
                if 'achiral' in chiral_types or 'aCh' in chiral_types:
                    if 'armchair' not in chiral_types:
                        chiral_types.append('armchair')
                    if 'zigzag' not in chiral_types:
                        chiral_types.append('zigzag')
                Ch_list = [Ch for Ch in Ch_list
                           if get_Ch_type(Ch) in chiral_types]
        return Ch_list


def generate_Ch_property_grid(compute=str, imax=10, **kwargs):
    """Generate a 2-dimensional,
    :math:`i_{\\mathrm{max}}\\times i_{\\mathrm{max}}` grid of
    nanotube properties.

    The property grid is indexed by :math:`(n, m)` chiralities
    for :math:`0\\le n\\le i_{\\mathrm{max}}` and
    :math:`0\\le m\\le i_{\\mathrm{max}}`.

    Parameters
    ----------
    compute_method : str
    imax : int

    Returns
    -------
    grid : ndarray

    """
    try:
        compute_func = \
            getattr(importlib.import_module('sknano.structures'), compute)
        grid = np.zeros((imax + 1, imax + 1)) - 1
        for n in range(imax + 1):
            for m in range(imax + 1):
                grid[n, m] = compute_func(n, m, **kwargs)
        return grid
    except AttributeError as e:
        print(e)
        return None


def get_chiral_type(Ch):
    """Identify the type of nanotube based on its chirality

    Parameters
    ----------
    Ch : {str, tuple}

    Returns
    -------
    str
        `armchair` for armchair tubes,
        `zigzag` for zigzag tubes,
        `chiral` for chiral tubes,

    """
    if isinstance(Ch, tuple):
        n, m = Ch
    else:
        n, m = get_chiral_indices_from_str(Ch)
    if n == m:
        return 'armchair'
    elif n != m and (n == 0 or m == 0):
        return 'zigzag'
    else:
        return 'chiral'

get_Ch_type = get_chiral_type


def map_Ch(Ch, compute=None, **kwargs):
    """Map compute function using `Ch` as input.

    Parameters
    ----------
    Ch : {str, tuple}
    compute : str

    Returns
    -------
    value of compute function.
    """
    try:
        compute_func = \
            getattr(importlib.import_module('sknano.structures'), compute)
        if isinstance(Ch, tuple):
            n, m = Ch
        else:
            n, m = get_chiral_indices_from_str(Ch)

        return compute_func(n, m, **kwargs)
    except AttributeError as e:
        print(e)
        return None


def get_Ch_map_from_data(Chlist, Chdata):
    pass


def get_chiral_indices(*args, check_type=True, **kwargs):
    """Parse the chiral indices `n` and `m` from a `vararg` `*args`, \
        which may be a `tuple` or 2 ints or from varkwargs `**kwargs`.

    """
    n = m = None
    try:
        n, m = args
    except ValueError:
        try:
            n, m = args[0]
        except IndexError:
            try:
                n, m = kwargs['Ch']
                del kwargs['Ch']
            except KeyError:
                n = kwargs['n']
                del kwargs['n']
                m = kwargs['m']
                del kwargs['m']
    if check_type:
        type_check_chiral_indices((n, m))
    return n, m, kwargs

get_Ch_indices = get_chiral_indices


def get_chiral_indices_from_str(Ch):
    """Return the chiral indicies `n` and `m`.

    Parameters
    ----------
    Ch : str
        string of the form 'NNMM' or "(n, m)". Extra spaces are acceptable.

    Returns
    -------
    sequence
        2-tuple of chiral indices `n` and `m`

    """
    try:
        return int(Ch[:2]), int(Ch[2:])
    except ValueError:
        try:
            return tuple([int(c) for c in re.split('[\(\),\s*]+', Ch)[1:-1]])
        except Exception as e:
            print(e)
            return None


def type_check_chiral_indices(*Ch):
    """Check type of chiral indices.

    Parameters
    ----------
    Ch : :class:`~python:tuple`

    Raises
    ------
    TypeError

    """
    n = m = None
    try:
        n, m = Ch
    except ValueError:
        try:
            n, m = Ch[0]
        except IndexError:
            raise ValueError('Expected an (n, m) tuple of chiral indices')

    if not (isinstance(n, numbers.Real) or n >= 0):
        raise TypeError('Expected an integer')
    if not (isinstance(m, numbers.Real) or m >= 0):
        raise TypeError('Expected an integer')
