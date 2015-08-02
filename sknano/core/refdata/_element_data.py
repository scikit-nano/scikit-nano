# -*- coding: utf-8 -*-
"""
=============================================================================
Reference element data (:mod:`sknano.core.refdata._element_data`)
=============================================================================

.. currentmodule:: sknano.core.refdata._element_data

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from sknano.core import loadobj
# from pkg_resources import resource_filename
import os

__all__ = ['element_data']


element_data = loadobj(os.path.join(os.path.dirname(__file__),
                                    'element_data.yaml'))

# dVDW = 3.35  # angstroms
