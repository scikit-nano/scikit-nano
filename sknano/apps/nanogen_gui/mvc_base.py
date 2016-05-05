# -*- coding: utf-8 -*-
"""
===============================================================================
NanoGen MVC base classes (:mod:`sknano.apps.nanogen_gui.mvc_base`)
===============================================================================

.. currentmodule:: sknano.apps.nanogen_gui.mvc_base

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'


from sknano.core.refdata import aCC, element_data

from .mvc_mixins import ObserverModelMixin, ViewControllerMixin

_r_CC_vdw = element_data['C']['VanDerWaalsRadius']

__all__ = ['GeneratorViewController', 'NanoStructureModelBase']


class GeneratorViewController(ViewControllerMixin):
    """Mixin Generator View Controller class."""
    def __init__(self, model=None, view=None, **kwargs):
        self.model = model
        self.view = view(self, self.model, **kwargs)
        self.model.notify_observers()
        self.view.show()

    def get_generator_parameters(self):
        return self.view.get_generator_parameters()


class NanoStructureModelBase(ObserverModelMixin):

    def __init__(self):
        super().__init__()
        self._observers = []
        self._bond = aCC
        self._element1 = self._element2 = 'C'

    @property
    def basis(self):
        return [self.element1, self.element2]

    @property
    def bond(self):
        return self._bond

    @bond.setter
    def bond(self, value):
        self._bond = value
        self.notify_observers()

    @property
    def element1(self):
        return self._element1

    @element1.setter
    def element1(self, value):
        self._element1 = value
        self.notify_observers()

    @property
    def element2(self):
        return self._element2

    @element2.setter
    def element2(self, value):
        self._element2 = value
        self.notify_observers()
