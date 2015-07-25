# -*- coding: utf-8 -*-
"""
=============================================
I/O functions (:mod:`sknano.core._io`)
=============================================

.. currentmodule:: sknano.core._io

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import json
import os
import re
import sys

try:
    import yaml
    try:
        from yaml import CLoader as Loader, CDumper as Dumper
    except ImportError:
        from yaml import Loader, Dumper
except ImportError:
    yaml = None


__all__ = ['get_fname', 'get_fpath', 'listdir_dirnames', 'listdir_fnames',
           'listdir', 'loadobj', 'dumpobj']


def get_fname(fname=None, ext=None, outpath=os.getcwd(), overwrite=False,
              add_fnum=True, fnum=None, verbose=False, **kwargs):
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
        If `fname.endswith(ext)` is `True` from the start, then `ext` will
        not be duplicately appended.
    outpath : str, optional
        Absolute or relative path for generated output file.
        Default is the absolute path returned by `os.getcwd()`.
    overwrite : bool, optional
        If `True`, overwrite an existing file if it has the same generated
        file path.
    add_fnum : bool, optional
        Append integer number to output file name, starting with **1**.
    fnum : {None, int}, optional
        Starting file number to append if `add_fnum` is `True`.

        .. note::

        If the generated file path exists and `overwrite` is False,
        setting this parameter has no effect.

    verbose : bool, optional
        Show verbose output.

    Returns
    -------
    fname : str
        Updated `fname`.

    """
    return get_fpath(fname=fname, ext=ext, outpath=outpath,
                     overwrite=overwrite, add_fnum=add_fnum, fnum=fnum,
                     fname_only=True, verbose=verbose)


def get_fpath(fname=None, ext=None, outpath=os.getcwd(), overwrite=False,
              add_fnum=True, fnum=None, include_fname=False, fname_only=False,
              verbose=False):
    """Generate absolute path to modified `fname`.

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
        If `fname.endswith(ext)` is `True` from the start, then `ext` will
        not be duplicately appended.
    outpath : str, optional
        Absolute or relative path for generated output file.
        Default is the absolute path returned by `os.getcwd()`.
    overwrite : bool, optional
        If `True`, overwrite an existing file if it has the same generated
        file path.
    add_fnum : bool, optional
        Append integer number to output file name, starting with **1**.
    fnum : {None, int}, optional
        Starting file number to append if `add_fnum` is `True`.

        .. note::

        If the generated file path exists and `overwrite` is False,
        setting this parameter has no effect.

    include_fname : bool, optional
        If `True`, return `(fpath, fname)` tuple.
    fname_only : bool, optional
        If `True`, return only `fname`.
    verbose : bool, optional
        Show verbose output.

    Returns
    -------
    fpath : str
        The concatenation of `outpath` followed by the updated `fname`.
    (fpath, fname) : tuple (only if `include_fname` is `True`)
        2-tuple of strings `(fpath, fname)`.
    fname : str (only if `fname_only` is `True`)
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

    if outpath is None:
        outpath = os.getcwd()

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
                        print('overwriting existing file: {}'.format(fname))
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
                    print('file exists: {}\n'.format(fpath))
                    print('Set `add_fnum=True` to generate unique\n' +
                          '`fname` or `overwrite=True` to overwrite\n' +
                          'existing file.')
                    fpath = None

        if verbose:
            print('Generated file name: {}'.format(fname))
            print('File path: {}'.format(fpath))

        if fname_only:
            return fname
        elif include_fname:
            return fpath, fname
        else:
            return fpath


def listdir_dirnames(path='.', filterfunc=None):
    """Return list of names of directories in the directory given by `path`.

    Parameters
    ----------
    path : :class:`~python:str`, optional
    filterfunc : `function`, optional

    Returns
    -------
    fnames : :class:`~python:list`
        :class:`~python:list` of names of directories in `path`.

    """
    return listdir(path, filterfunc=filterfunc,
                   filter_dirnames=filterfunc is not None)[0]


def listdir_fnames(path='.', filterfunc=None):
    """Return list of names of files in the directory given by `path`.

    Parameters
    ----------
    path : :class:`~python:str`, optional
    filterfunc : `function`, optional

    Returns
    -------
    fnames : :class:`~python:list`
        :class:`~python:list` of names of files in `path`.

    """
    return listdir(path, filterfunc=filterfunc,
                   filter_fnames=filterfunc is not None)[-1]


def listdir(path='.', filterfunc=None, filter_dirnames=False,
            filter_fnames=False):
    """Return a tuple of the names of the directories and files in the
    directory given by `path`.

    Parameters
    ----------
    path : :class:`~python:str`, optional
    filterfunc : `function`, optional
    filter_dirnames : :class:`~python:bool`, optional
    filter_fnames : :class:`~python:bool`, optional

    Returns
    -------
    (dirnames, fnames) : :class:`~python:tuple`
        :class:`~python:tuple` of names of directories and files in `path`.

    """
    l = os.listdir(path)
    dirnames = [name for name in l if os.path.isdir(os.path.join(path, name))]
    fnames = [name for name in l if os.path.isfile(os.path.join(path, name))]
    if filter_dirnames and filterfunc is not None:
        dirnames = list(filter(filterfunc, dirnames))
    if filter_fnames and filterfunc is not None:
        fnames = list(filter(filterfunc, fnames))
    return dirnames, fnames


def loadobj(fn, *args, **kwargs):
    with open(fn) as fp:
        if fn.lower().endswith(("yaml", "yml")):
            if "Loader" not in kwargs:
                kwargs["Loader"] = Loader
            return yaml.load(fp, *args, **kwargs)
        else:
            return json.load(fp, *args, **kwargs)


def dumpobj(obj, fn, *args, **kwargs):
    with open(fn, 'wt') as fp:
        if fn.lower().endswith(("yaml", "yml")):
            if "Dumper" not in kwargs:
                kwargs["Dumper"] = Dumper
            yaml.dump(obj, fp, *args, **kwargs)
        else:
            json.dump(obj, fp, *args, **kwargs)
