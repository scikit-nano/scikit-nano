# -*- coding: utf-8 -*-
"""
==============================================================================
Crystallography helper functions (:mod:`sknano.core.crystallography._extras`)
==============================================================================

.. currentmodule:: sknano.core.crystallography._extras

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import inspect
import numpy as np

__all__ = ['pymatgen_structure']


def pymatgen_structure(*args, classmethod=None, scaling_matrix=None, **kwargs):
    """Helper function to generate a `pymatgen` Structure object.

    Parameters
    ----------
    *args : variable number of positional arguments
    classmethod : {None, str}, optional
    scaling_matrix : {class:`~python:int` or :class:`~python:list`}
    **kwargs : variable number of keyword arguments

    Returns
    -------
    :class:`pymatgen:Structure`

    """
    try:
        from pymatgen import Structure
    except ImportError as e:
        print(e)
    else:
        constructor = None
        if classmethod is None:
            constructor = Structure
        else:
            constructor = getattr(Structure, classmethod, None)

        pmg_sig = inspect.signature(constructor)
        bound_sig = pmg_sig.bind(*args, **kwargs)
        structure = constructor(*bound_sig.args, **bound_sig.kwargs)
        if scaling_matrix is not None and \
                isinstance(scaling_matrix,
                           (int, float, tuple, list, np.ndarray)):
            structure.make_supercell(scaling_matrix)
        return structure
