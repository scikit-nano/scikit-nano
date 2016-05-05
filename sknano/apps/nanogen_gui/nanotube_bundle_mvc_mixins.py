# -*- coding: utf-8 -*-
"""
=======================================================================================
Nanotube Bundle MVC mixins (:mod:`sknano.apps.nanogen_gui.nanotube_bundle_mvc_mixins`)
=======================================================================================

.. currentmodule:: sknano.apps.nanogen_gui.nanotube_bundle_mvc_mixins

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'


try:
    from PyQt5.QtCore import pyqtSlot
except ImportError:
    from PyQt4.QtCore import pyqtSlot


__all__ = ['NanotubeBundleModelMixin', 'NanotubeBundleViewMixin']


class NanotubeBundleModelMixin:
    @property
    def l1(self):
        return self.structure.lattice.a1.length

    @property
    def l2(self):
        return self.structure.lattice.a2.length

    @property
    def n1(self):
        return self.structure.n1

    @n1.setter
    def n1(self, value):
        self.structure.n1 = value
        self.notify_observers()

    @property
    def n2(self):
        return self.structure.n2

    @n2.setter
    def n2(self, value):
        self.structure.n2 = value
        self.notify_observers()


class NanotubeBundleViewMixin:

    @pyqtSlot(int)
    def on_bundle_generator_check_box_stateChanged(self, value):
        [spin_box.setReadOnly(False if value else True)
         for spin_box in (self.bundle_n1_spin_box, self.bundle_n2_spin_box)]
        self.update_app_view()

    @pyqtSlot(int)
    def on_bundle_n1_spin_box_valueChanged(self, value):
        self.model.n1 = value

    @pyqtSlot(int)
    def on_bundle_n2_spin_box_valueChanged(self, value):
        self.model.n2 = value
