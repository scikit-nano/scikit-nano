# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

from .collections import *
from .io import *
from .itertools import *
from .meta import *
from .strings import *

__all__ = ['components', 'dimensions', 'xyz', 'xyz_axes']

components = dimensions = xyz = xyz_axes = ('x', 'y', 'z')

__all__ += [s for s in dir() if not s.startswith('_')]
