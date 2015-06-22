# -*- coding: utf-8 -*-
"""
==============================================================================
View mixin classes (:mod:`sknano.apps.nanogen_gui._view_mixins`)
==============================================================================

.. currentmodule:: sknano.apps.nanogen_gui._view_mixins

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

try:
    # from PyQt4.QtCore import pyqtSlot
    from PyQt5.QtCore import pyqtSlot
except ImportError as e:
    print(e)

__all__ = ['MainWindowViewMixin', 'SWNTViewMixin', 'MWNTViewMixin',
           'GrapheneViewMixin', 'FullereneViewMixin']


class MainWindowViewMixin:
    """Mixin class for main window."""

    @pyqtSlot(str)
    def on_element1_comboBox_currentIndexChanged(self, value):
        self.model.element1 = str(value)

    @pyqtSlot(str)
    def on_element2_comboBox_currentIndexChanged(self, value):
        self.model.element2 = str(value)

    @pyqtSlot()
    def on_bond_doubleSpinBox_editingFinished(self):
        self.model.bond = self.bond_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_bond_doubleSpinBox_valueChanged(self, value):
        self.model.bond = value


class SWNTViewMixin:
    """Mixin class for nanotube classes."""
    @pyqtSlot()
    def on_swnt_n_spinBox_editingFinished(self):
        self.model.n = self.swnt_n_spinBox.value()

    @pyqtSlot(int)
    def on_swnt_n_spinBox_valueChanged(self, value):
        self.model.n = value

    @pyqtSlot()
    def on_swnt_m_spinBox_editingFinished(self):
        self.model.m = self.swnt_m_spinBox.value()

    @pyqtSlot(int)
    def on_swnt_m_spinBox_valueChanged(self, value):
        self.model.m = value

    @pyqtSlot()
    def on_swnt_nz_doubleSpinBox_editingFinished(self):
        self.model.nz = self.swnt_nz_doubleSpinBox.value()

    @pyqtSlot(int)
    def on_swnt_nz_doubleSpinBox_valueChanged(self, value):
        self.model.nz = value

    @pyqtSlot()
    def on_swnt_Lz_doubleSpinBox_editingFinished(self):
        self.model.Lz = self.swnt_Lz_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_swnt_Lz_doubleSpinBox_valueChanged(self, value):
        self.model.Lz = value

    @pyqtSlot()
    def on_swnt_bundle_nx_spinBox_editingFinished(self):
        self.model.nx = self.nx_spinBox.value()

    @pyqtSlot(int)
    def on_swnt_bundle_nx_spinBox_valueChanged(self, value):
        self.model.nx = value

    @pyqtSlot()
    def on_swnt_bundle_ny_spinBox_editingFinished(self):
        self.model.ny = self.ny_spinBox.value()

    @pyqtSlot(int)
    def on_swnt_bundle_ny_spinBox_valueChanged(self, value):
        self.model.ny = value

    # @pyqtSlot(int)
    # def on_swnt_tab_generator_buttonGroup_buttonClicked(self):
    #     if self.unrolled_swnt_generator_radioButton.isChecked():
    #         if self.bundle_generator_checkBox.isChecked():
    #             self.bundle_generator_checkBox.setCheckState(False)
    #         self.bundle_generator_checkBox.setCheckable(False)
    #     else:
    #         self.bundle_generator_checkBox.setCheckable(True)


class MWNTViewMixin:

    @pyqtSlot()
    def on_mwnt_bundle_nx_spinBox_editingFinished(self):
        self.model.nx = self.mwnt_bundle_nx_spinBox.value()

    @pyqtSlot(int)
    def on_mwnt_bundle_nx_spinBox_valueChanged(self, value):
        self.model.nx = value

    @pyqtSlot()
    def on_mwnt_bundle_ny_spinBox_editingFinished(self):
        self.model.ny = self.mwnt_bundle_ny_spinBox.value()

    @pyqtSlot(int)
    def on_mwnt_bundle_ny_spinBox_valueChanged(self, value):
        self.model.ny = value

    @pyqtSlot()
    def on_mwnt_Lz_spinBox_editingFinished(self):
        self.model.Lz = self.mwnt_Lz_doubleSpinBox.value()

    @pyqtSlot(int)
    def on_mwnt_Lz_doubleSpinBox_valueChanged(self, value):
        self.model.Lz = value


class GrapheneViewMixin:
    """Mixin class for graphene."""

    # @pyqtSlot(int)
    # def on_stacking_order_buttonGroup_buttonClicked(self):
    #     if self.nlayer_spinBox.value() == 1:
    #         pass

    @pyqtSlot()
    def on_armchair_edge_length_doubleSpinBox_editingFinished(self):
        self.model.armchair_edge_length = \
            self.armchair_edge_length_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_armchair_edge_length_doubleSpinBox_valueChanged(self, value):
        self.model.armchair_edge_length = value

    @pyqtSlot()
    def on_zigzag_edge_length_doubleSpinBox_editingFinished(self):
        self.model.zigzag_edge_length = \
            self.zigzag_edge_length_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_zigzag_edge_length_doubleSpinBox_valueChanged(self, value):
        self.model.zigzag_edge_length = value

    @pyqtSlot()
    def on_nlayers_spinBox_editingFinished(self):
        self.model.nlayers = self.nlayers_spinBox.value()
        #self.update_stacking_order_buttonGroup()

    @pyqtSlot(int)
    def on_nlayers_spinBox_valueChanged(self, value):
        self.model.nlayers = value
        # self.update_stacking_order_buttonGroup()

    # def update_stacking_order_buttonGroup(self):
    #     if self.nlayers_spinBox.value() == 1:
    #         for rb in (self.AA_stacking_radioButton,
    #                    self.AB_stacking_radioButton):
    #             rb.setCheckable(False)
    #     else:
    #         for rb in (self.AA_stacking_radioButton,
    #                    self.AB_stacking_radioButton):
    #             rb.setCheckable(True)


class FullereneViewMixin:
    pass
