# -*- coding: utf-8 -*-
"""
=================================================================
nanogen helper functions (:mod:`sknano.nanogen._nanogen_funcs`)
=================================================================

.. currentmodule:: sknano.nanogen._nanogen_funcs

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numpy as np

from ._nanotubes import Nanotube
from ..tools import chiral_type_name_mappings as Ch_types

__all__ = ['generate_Ch_list', 'generate_Ch_property_grid']


def generate_Ch_list(ns=None, ni=None, nf=None, dn=None,
                     ms=None, mi=None, mf=None, dm=None,
                     chiral_type=None, handedness=None,
                     echo_zsh_str=False):
    """Generate a list of :math:`(n, m)` chiralities.

    Parameters
    ----------
    ns, ms : sequence, optional
        list of :math:`n` and :math:`m` chiral indices
    ni, nf, dn : int, optional
        :math:`(n_i, n_f, \\Delta n)` denote the `start`, `stop`, and
        `step` parameters passed to the numpy function `np.arange`,
        to generate an array of evenly spaced :math:`n` chiral indices.
    mi, mf, dm : int
        :math:`(m_i, m_f, \\Delta m)` denote the `start`, `stop`, and
        `step` parameters passed to the numpy function `np.arange`,
        to generate an array of evenly spaced :math:`m` chiral indices.
    chiral_type : {'all', 'achiral', 'chiral', 'armchair', 'zigzag'}, optional
    handedness : {'all', 'left', 'right'}, optional
    echo_zsh_str : bool, optional

    Returns
    -------
    Ch_list : sequence
        list of chiralities

    """
    Ch_list = []
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
            if chiral_type in (Ch_types.keys() + Ch_types.values()):
                if chiral_type in ('achiral', 'aCh'):
                    for n in ns:
                        Ch_list.append((n, n))
                    for n in ns:
                        Ch_list.append((n, 0))
                elif chiral_type in ('armchair', 'AC'):
                    for n in ns:
                        Ch_list.append((n, n))
                elif chiral_type in ('zigzag', 'ZZ'):
                    for n in ns:
                        Ch_list.append((n, 0))
                else:
                    for n in ns:
                        m = 1
                        while m < n:
                            if handedness == 'left':
                                Ch_list.append((m, n))
                            elif handedness == 'right':
                                Ch_list.append((n, m))
                            else:
                                Ch_list.append((n, m))
                                Ch_list.append((m, n))
                            m += 1
            else:
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
                        if chiral_type not in \
                                ('chiral', 'Ch', 'armchair', 'AC'):
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
                    elif (n == m):
                        Ch = (n, n)
                        if chiral_type not in \
                                ('chiral', 'Ch', 'zigzag', 'ZZ') \
                                and Ch not in Ch_list:
                            Ch_list.append(Ch)
                    else:
                        nm = (n, m)
                        mn = (m, n)
                        if chiral_type not in ('achiral', 'aCh', 'armchair',
                                               'AC', 'zigzag', 'ZZ'):
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

        return Ch_list


def generate_Ch_property_grid(compute_method=str, imax=10, **kwargs):
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
        nanotube_compute = getattr(Nanotube, compute_method)
        grid = np.zeros((imax + 1, imax + 1)) - 1
        for ni in xrange(imax + 1):
            for mi in xrange(imax + 1):
                grid[ni, mi] = nanotube_compute(n=ni, m=mi, **kwargs)
        return grid
    except AttributeError as e:
        print(e)
        return None
