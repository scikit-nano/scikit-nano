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

__all__ = ['NanotubeViewMixin', 'GrapheneViewMixin']


class NanotubeViewMixin:
    """Mixin class for nanotube classes."""
    @pyqtSlot()
    def on_n_spinBox_editingFinished(self):
        self.model.n = self.n_spinBox.value()

    @pyqtSlot(int)
    def on_n_spinBox_valueChanged(self, value):
        self.model.n = value

    @pyqtSlot()
    def on_m_spinBox_editingFinished(self):
        self.model.m = self.m_spinBox.value()

    @pyqtSlot(int)
    def on_m_spinBox_valueChanged(self, value):
        self.model.m = value

    @pyqtSlot()
    def on_nanotube_bond_doubleSpinBox_editingFinished(self):
        self.model.nanotube_bond = self.nanotube_bond_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_nanotube_bond_doubleSpinBox_valueChanged(self, value):
        self.model.nanotube_bond = value

    @pyqtSlot()
    def on_nx_spinBox_editingFinished(self):
        self.model.nx = self.nx_spinBox.value()

    @pyqtSlot(int)
    def on_nx_spinBox_valueChanged(self, value):
        self.model.nx = value

    @pyqtSlot()
    def on_ny_spinBox_editingFinished(self):
        self.model.ny = self.ny_spinBox.value()

    @pyqtSlot(int)
    def on_ny_spinBox_valueChanged(self, value):
        self.model.ny = value

    @pyqtSlot()
    def on_nz_spinBox_editingFinished(self):
        self.model.nz = self.nz_spinBox.value()

    @pyqtSlot(int)
    def on_nz_spinBox_valueChanged(self, value):
        self.model.nz = value

    @pyqtSlot(int)
    def on_nanotube_generator_buttonGroup_buttonClicked(self):
        if self.unrolled_nanotube_generator_radioButton.isChecked():
            if self.bundle_generator_checkBox.isChecked():
                self.bundle_generator_checkBox.setCheckState(False)
            self.bundle_generator_checkBox.setCheckable(False)
        else:
            self.bundle_generator_checkBox.setCheckable(True)

    @pyqtSlot(str)
    def on_nanotube_element1_comboBox_currentIndexChanged(self, value):
        self.model.nanotube_element1 = str(value)

    @pyqtSlot(str)
    def on_nanotube_element2_comboBox_currentIndexChanged(self, value):
        self.model.nanotube_element2 = str(value)

    @pyqtSlot(str)
    def on_graphene_element1_comboBox_currentIndexChanged(self, value):
        self.model.graphene_element1 = str(value)

    @pyqtSlot(str)
    def on_graphene_element2_comboBox_currentIndexChanged(self, value):
        self.model.graphene_element2 = str(value)

    @pyqtSlot()
    def on_Lz_doubleSpinBox_editingFinished(self):
        self.model.Lz = self.Lz_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_Lz_doubleSpinBox_valueChanged(self, value):
        self.model.Lz = value


class GrapheneViewMixin(object):
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

    @pyqtSlot()
    def on_graphene_bond_doubleSpinBox_editingFinished(self):
        self.model.graphene_bond = self.graphene_bond_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_graphene_bond_doubleSpinBox_valueChanged(self, value):
        self.model.graphene_bond = value
