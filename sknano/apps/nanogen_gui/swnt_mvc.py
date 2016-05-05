# -*- coding: utf-8 -*-
"""
===============================================================================
SWNT MVC classes (:mod:`sknano.apps.nanogen_gui.swnt_mvc`)
===============================================================================

.. currentmodule:: sknano.apps.nanogen_gui.swnt_mvc

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

try:
    from PyQt5.QtCore import pyqtSlot
    from PyQt5.QtWidgets import QMainWindow
    from ._pyqt5_ui_swnt_generator import Ui_SWNTGenerator
except ImportError:
    from PyQt4.QtCore import pyqtSlot
    from PyQt4.QtGui import QMainWindow
    from ._pyqt4_ui_swnt_generator import Ui_SWNTGenerator

from sknano.core.structures import SWNT, compute_L, compute_T

from .mvc_base import GeneratorViewController, NanoStructureModelBase
from .mvc_mixins import NanoStructureViewMixin
from .nanotube_bundle_mvc_mixins import NanotubeBundleModelMixin, \
    NanotubeBundleViewMixin

__all__ = ['SWNTModelMixin', 'SWNTModel', 'SWNTViewMixin', 'SWNTGeneratorView',
           'SWNTGeneratorController']


class SWNTModelMixin:
    @property
    def Ch(self):
        """Chiral indices :math:`(n, m)`"""
        return self.n, self.m

    @property
    def n(self):
        """Chiral index :math:`n`"""
        return self.structure.n

    @n.setter
    def n(self, value):
        self.structure.n = value
        self.notify_observers()

    @property
    def m(self):
        """Chiral index :math:`m`"""
        return self.structure.m

    @m.setter
    def m(self, value):
        self.structure.m = value
        self.notify_observers()

    @property
    def L(self):
        return compute_L(self.Ch, n3=self.n3, bond=self.bond)

    @L.setter
    def L(self, value):
        n3 = value / compute_T(self.Ch, bond=self.bond, length=True)
        if not self.structure.fix_L:
            n3 = int(n3)
        self.structure.n3 = n3
        self.notify_observers()

    @property
    def n3(self):
        return self.structure.n3

    @n3.setter
    def n3(self, value):
        self.structure.n3 = value
        self.notify_observers()

    @property
    def fix_L(self):
        return self.structure.fix_L

    @fix_L.setter
    def fix_L(self, value):
        self.structure.fix_L = value
        self.notify_observers()


class SWNTViewMixin:
    """Mixin class for nanotube classes."""
    @pyqtSlot(int)
    def on_swnt_n_spin_box_valueChanged(self, value):
        self.model.n = value

    @pyqtSlot(int)
    def on_swnt_m_spin_box_valueChanged(self, value):
        self.model.m = value

    @pyqtSlot(float)
    def on_swnt_n3_double_spin_box_valueChanged(self, value):
        self.model.n3 = value

    @pyqtSlot(float)
    def on_swnt_L_double_spin_box_valueChanged(self, value):
        self.model.L = value

    @pyqtSlot(int)
    def on_swnt_fix_L_check_box_stateChanged(self, value):
        if value:
            self.model.fix_L = True
            self.swnt_n3_double_spin_box.setReadOnly(True)
            self.swnt_L_double_spin_box.setReadOnly(False)
        else:
            self.model.fix_L = False
            self.swnt_n3_double_spin_box.setReadOnly(False)
            self.swnt_L_double_spin_box.setReadOnly(True)
        self.update_app_view()


class SWNTModel(NanoStructureModelBase, SWNTModelMixin,
                NanotubeBundleModelMixin):
    def __init__(self):
        super().__init__()
        self.structure = SWNT((10, 10), basis=self.basis, bond=self.bond,
                              n1=1, n2=1, n3=1, bundle_packing='hcp')
        self.notify_observers()


class SWNTGeneratorView(QMainWindow, Ui_SWNTGenerator, NanoStructureViewMixin,
                        SWNTViewMixin, NanotubeBundleViewMixin):

    def __init__(self, controller=None, model=None, parent=None):
        self.controller = controller
        self.model = model
        self.parent = parent
        model.register_observer(self)
        super().__init__(parent=parent)
        self.setupUi(self)

    def get_generator_parameters(self):
        kwargs = {}
        element1 = str(self.element1_combo_box.itemText(
                       self.element1_combo_box.currentIndex()))
        element2 = str(self.element2_combo_box.itemText(
                       self.element2_combo_box.currentIndex()))
        kwargs['basis'] = [element1, element2]
        kwargs['bond'] = self.bond_double_spin_box.value()

        kwargs['n'] = self.swnt_n_spin_box.value()
        kwargs['m'] = self.swnt_m_spin_box.value()
        kwargs['fix_L'] = self.swnt_fix_L_check_box.isChecked()
        kwargs['L'] = self.swnt_L_double_spin_box.value()
        kwargs['n3'] = self.swnt_n3_double_spin_box.value()

        kwargs['generator_class'] = 'SWNTGenerator'
        if self.bundle_generator_check_box.isChecked():
            kwargs['bundle_packing'] = 'hcp'
            kwargs['n1'] = self.bundle_n1_spin_box.value()
            kwargs['n2'] = self.bundle_n2_spin_box.value()
            # kwargs['Ntubes'] = self.model.structure.Ntubes

        return kwargs

    def update_app_view(self):
        self.elements_bond_label.setText('-'.join((self.model.element1,
                                                   self.model.element2,
                                                   ' bond =')))
        self.bond_double_spin_box.setValue(self.model.bond)

        self.swnt_n_spin_box.setValue(self.model.n)
        self.swnt_m_spin_box.setValue(self.model.m)

        self.swnt_n3_double_spin_box.blockSignals(True)
        self.swnt_L_double_spin_box.blockSignals(True)
        self.swnt_L_double_spin_box.setValue(self.model.L)
        self.swnt_n3_double_spin_box.setValue(self.model.n3)
        self.swnt_n3_double_spin_box.blockSignals(False)
        self.swnt_L_double_spin_box.blockSignals(False)

        self.bundle_n1_spin_box.setValue(self.model.n1)
        self.bundle_n2_spin_box.setValue(self.model.n2)

        # self.bundle_l1_line_edit.setText('{:.4f} Å'.format(self.model.l1))
        # self.bundle_l2_line_edit.setText('{:.4f} Å'.format(self.model.l2))
        self.parent.update_app_view()


class SWNTGeneratorController(GeneratorViewController):
    def __init__(self, model=None, **kwargs):
        super().__init__(model=model, view=SWNTGeneratorView, **kwargs)
