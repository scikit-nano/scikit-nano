# -*- coding: utf-8 -*-
"""
==============================================
Helper functions (:mod:`sknano.tools._funcs`)
==============================================

.. currentmodule:: sknano.tools._funcs

"""
from __future__ import division, print_function, absolute_import

from fractions import gcd
import os
import re
import sys

import numpy as np

from ._luts import chiral_type_name_mappings as Ch_types

__all__ = ['cmp_Ch', 'filter_Ch', 'filter_Ch_list', 'filter_key_type_mappings',
           'generate_Ch_list', 'get_Ch_indices', 'get_Ch_type',
           'get_fpath', 'plural_word_check', 'rotation_matrix', 'totient_func']

filter_key_type_mappings = {}
filter_key_type_mappings['Ch_type'] = str
for k in ('even', 'odd'):
    filter_key_type_mappings[k + '_only'] = bool

for k in ('min_index', 'max_index', 'min_n', 'max_n', 'min_m', 'max_m',
          'n', 'm'):
    filter_key_type_mappings[k] = int


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
    elif isinstance(Ch1, str) and isinstance(Ch2, str):
        n1, m1 = get_Ch_indices(Ch1)
        Ch1_type = get_Ch_type(Ch1)
        n2, m2 = get_Ch_indices(Ch2)
        Ch2_type = get_Ch_type(Ch2)

    if Ch1_type == 'AC' and Ch2_type == 'ZZ':
        return 1
    elif Ch1_type == 'ZZ' and Ch2_type == 'AC':
        return -1
    #if dt1 > dt2:
    #    return 1
    #elif dt1 < dt2:
    #    return -1
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
        n, m = get_Ch_indices(Ch)
        this_Ch_type = get_Ch_type(Ch)

    if even_only:
        if n % 2 == 0 and m % 2 == 0:
            return True
        else:
            return False

    elif odd_only:
        if this_Ch_type != 'ZZ':
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
        if Ch_type in Ch_types:
            Ch_type = Ch_types[Ch_type]

        if this_Ch_type == Ch_type:
            return True
        else:
            return False

    elif min_index is not None and isinstance(min_index, int):
        if this_Ch_type != 'ZZ':
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


def filter_Ch_list(Ch_list, **kwargs):
    """Filter for filtering list of chiralities against set of conditions.

    Parameters
    ----------
    Ch_list : sequence
        list of chiralities
    kwargs : dict, optional
        dictionary of conditions to test each chirality against.

    Returns
    -------
    sequence
        list of chiralities which passed conditions

    """
    return [Ch for Ch in Ch_list if filter_Ch(Ch, **kwargs)]


def get_Ch_indices(Ch):
    """Return the chiral indicies `n` and `m`.

    Parameters
    ----------
    Ch : str

    Returns
    -------
    sequence
        list of chiral indices `n` and `m`

    """
    return [int(c) for c in re.split('[\(\),\s*]+', Ch)[1:-1]]


#def get_nearest_neighbor_counts(adata=None, datafile=None,
#                                filter_elements=None, invert=False):
#    if adata is not None:
#        atom_data = adata
#    elif datafile is not None and isinstance(datafile, str):
#        atom_data = get_atom_data_asarray(datafile)
#
#    #atom_ids = get_atom_ids_asarray(adata=atom_data)
#
#    filtered_coords = \
#        get_filtered_coords(adata=atom_data, filter_elements=vacIDs)
#    NNtree = spatial.KDTree(filtered_coords)
#    print(NNtree.data)
#
#    # compute things like the distance for which 50% of vacancies
#    # lie within a distance r of each other, etc.
#
#    #return None
#    pass

#def get_nearest_neighbor_comp_array_dict(vacIDs):
#    nearest_neighbor_comp_dict = OrderedDict()
#    vac_coords = \
#         get_filtered_coords(adata=self.atom_data, filter_elements=vacIDs)
#    #vac_coords = self.get_filtered_coords(vacIDs)
#    for comp in self.position_components:
#        vac_distances = self.get_vac_distances(vac_coords, comp)
#        nearest_neighbor_array = vac_distances[np.where(vac_distances == \
#                                                        vac_distances.min())]
#        nearest_neighbor_comp_dict[comp] = nearest_neighbor_array
#
#    return nearest_neighbor_comp_dict


def get_Ch_type(Ch):
    """Identify the type of nanotube based on its chirality

    Parameters
    ----------
    Ch : {str, tuple}

    Returns
    -------
    str
        `AC` for armchair tubes,
        `ZZ` for zigzag tubes,
        `Ch` for chiral tubes,

    """
    if isinstance(Ch, tuple):
        n, m = Ch
    else:
        n, m = get_Ch_indices(Ch)
    if n == m:
        return 'AC'
    elif n != m and (n == 0 or m == 0):
        return 'ZZ'
    else:
        return 'Ch'


def generate_Ch_list(ns=None, ni=None, nf=None, dn=None,
                     ms=None, mi=None, mf=None, dm=None,
                     chiral_type=None, handedness=None,
                     echo_zsh_str=False):
    """Generate a list of chiralities.

    Parameters
    ----------
    ns, ms : sequence
        list of n and m chiral indices
    ni, nf, dn : int
    mi, mf, dm : int
    chiral_type : {None, 'all', 'achiral', 'chiral', 'armchair', 'zigzag'}
    handedness : {None, 'all', 'left', 'right'}
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


def get_fpath(fname=None, ext=None, outpath=os.getcwd(), overwrite=False,
              add_fnum=True, fnum=None, include_fname=False, fname_only=False,
              verbose=False):
    """Generate modified `fname` string based on chosen parameters.

    Parameters
    ----------
    fname : str
        Name of file, with or without an extension.
    ext : str, optional
        File extension to append to `fname`. If `ext` is None,
        then `fname` is analyzed to see if it likely already has an
        extension. An extension is set to the
        last element in the list of strings returned by
        `fname.split('.')` **if** this list has more than 1 element.
        Otherwise, `ext` will be set to an empty string `''`.
        If `ext` is not None and is a valid string,
        then `fname` is analyzed to see if it already ends with `ext`.
        If `fname.endswith(ext)` is True from the start, then `ext` will
        not be duplicately appended.
    outpath : str, optional
        Absolute or relative path for generated output file.
        Default is the absolute path returned by `os.getcwd()`.
    overwrite : bool, optional
        If True, overwrite an existing file if it has the same generated
        file path.
    add_fnum : bool, optional
        Append integer number to output file name, starting with **1**.
    fnum : {None, int}, optional
        Starting file number to append if `add_fnum` is True.

        .. note::

        If the generated file path exists and `overwrite` is False,
        setting this parameter has no effect.

    include_fname : bool, optional
        If True, return `(fpath, fname)` tuple.
    fname_only : bool, optional
        If True, return only `fname`.
    verbose : bool, optional
        Show verbose output.

    Returns
    -------
    fpath : str
        The concatenation of `outpath` followed by the updated `fname`.
    (fpath, fname) : tuple (only if `include_fname` is True)
        2-tuple of strings `(fpath, fname)`.
    fname : str (only if `fname_only` is True)
        Updated `fname`.

    """
    f = None
    if fname is None or fname == '':
        error_msg = '`fname` must be a string at least 1 character long.'
        if fname is None:
            raise TypeError(error_msg)
        else:
            raise ValueError(error_msg)
    else:
        f = fname
        fsplit = f.split('.')
        if ext is None:
            if len(fsplit) > 1:
                ext = '.' + fsplit[-1]
            else:
                ext = ''
        else:
            # check if extension already starts with a '.'
            if not ext.startswith('.'):
                ext = '.' + ext
            # check if file name already ends with extension.
            if f.split('.')[-1] != ext.split('.')[-1]:
                f += ext

    if add_fnum:
        fname = re.split(ext, f)[0]
        if fnum is not None:
            f = fname + '-{:d}'.format(fnum) + ext
        else:
            f = fname + '-1' + ext

    fpath = None

    try:
        os.makedirs(outpath)
    except OSError:
        if os.path.isdir(outpath):
            pass
        else:
            outpath = os.curdir
    finally:
        fname = f
        fpath = os.path.join(outpath, fname)
        if os.path.isfile(fpath):
            if overwrite:
                try:
                    os.remove(fpath)
                except OSError as e:
                    print(e)
                    sys.exit(1)
                else:
                    if verbose:
                        print(u'overwriting existing file: {}'.format(fname))
            else:
                if add_fnum:
                    while os.path.isfile(fpath):
                        fname = \
                            '-'.join(re.split('-', re.split(ext, f)[0])[:-1])
                        fnum = re.split('-', re.split(ext, f)[0])[-1]
                        f = fname + '-' + str(int(fnum) + 1) + ext
                        fpath = os.path.join(outpath, f)
                    fname = f
                else:
                    print(u'file exists: {}\n'.format(fpath) +
                          u'Set `add_fnum=True` to generate unique `fname`\n'
                          u'or `overwrite=True` to overwrite existing '
                          u'file.')
                    fpath = None

        if verbose:
            print(u'Generated file name: {}'.format(fname))
            print(u'File path: {}'.format(fpath))

        if fname_only:
            return fname
        elif include_fname:
            return fpath, fname
        else:
            return fpath


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


def rotation_matrix(angle, rot_axis='z', deg2rad=False):
    """Generate a rotation matrix.

    Parameters
    ----------
    angle : float
        rotation angle in radians
    rot_axis : {'x', 'y', 'z'}, optional
        rotation axis
    deg2rad : bool, optional
        Angle is in degrees and needs to be converted to radians

    Returns
    -------
    ndarray
        rotation matrix

    """
    if deg2rad:
        angle = np.radians(angle)
    if rot_axis == 'x':
        return np.array([[1, 0, 0],
                         [0, np.cos(angle), -np.sin(angle)],
                         [0, np.sin(angle), np.cos(angle)]])
    elif rot_axis == 'y':
        return np.array([[np.cos(angle), 0, np.sin(angle)],
                         [0, 1, 0],
                         [-np.sin(angle), 0, np.cos(angle)]])
    else:
        return np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0],
                         [0, 0, 1]])


def totient_func(n=int):
    """Compute the totatives of :math:`n`.

    Parameters
    ----------
    n : int

    Returns
    -------
    int : the number of totatives of :math:`n`

    """
    if not isinstance(n, int) and int(n) != n:
        raise ValueError('n must be an integer')
    n = int(n)

    tots = 0
    for k in xrange(1, n + 1):
        if gcd(n, k) == 1:
            tots += 1

    return tots
