# -*- coding: utf-8 -*-
"""
===============================================================================
NanoGen MVC mixins (:mod:`sknano.apps.nanogen_gui.mvc_mixins`)
===============================================================================

.. currentmodule:: sknano.apps.nanogen_gui.mvc_mixins

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

try:
    from PyQt5.QtCore import pyqtSlot
except ImportError:
    from PyQt4.QtCore import pyqtSlot

__all__ = ['ObserverModelMixin', 'ViewControllerMixin',
           'NanoStructureViewMixin']


class ObserverModelMixin:
    """Observer model mixin for registering/removing/notifying observers."""
    def register_observer(self, observer):
        self._observers.append(observer)

    def remove_observer(self, observer):
        self._observers.remove(observer)

    def notify_observers(self):
        for observer in self._observers[:]:
            observer.update_app_view()


class ViewControllerMixin:
    def refresh_view(self):
        """Refresh `NanoGenView`."""
        self.view.update_app_view()


class NanoStructureViewMixin:
    """Mixin class for main window."""
    @pyqtSlot(str)
    def on_element1_combo_box_currentIndexChanged(self, value):
        self.model.element1 = str(value)

    @pyqtSlot(str)
    def on_element2_combo_box_currentIndexChanged(self, value):
        self.model.element2 = str(value)

    @pyqtSlot(float)
    def on_bond_double_spin_box_valueChanged(self, value):
        self.model.bond = value
