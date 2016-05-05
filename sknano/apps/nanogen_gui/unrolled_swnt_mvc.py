# -*- coding: utf-8 -*-
"""
===============================================================================
Unrolled SWNT MVC classes (:mod:`sknano.apps.nanogen_gui.unrolled_swnt_mvc`)
===============================================================================

.. currentmodule:: sknano.apps.nanogen_gui.unrolled_swnt_mvc

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'


class UnrolledSWNTModelMixin:

    @property
    def n1(self):
        return self.structure.n1

    @n1.setter
    def n1(self, value):
        self.structure.n1 = value
        self.notify_observers()
