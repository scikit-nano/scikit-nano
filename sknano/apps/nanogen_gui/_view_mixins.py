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

__all__ = ['NanoGenViewMixin', 'SWNTViewMixin', 'MWNTViewMixin',
           'BundleViewMixin', 'GrapheneViewMixin', 'FullereneViewMixin',
           'BulkStructureViewMixin']


class NanoGenViewMixin:
    """Mixin class for main window."""
    @pyqtSlot(str)
    def on_element1_combo_box_currentIndexChanged(self, value):
        self.model.element1 = str(value)

    @pyqtSlot(str)
    def on_element2_combo_box_currentIndexChanged(self, value):
        self.model.element2 = str(value)

    # @pyqtSlot()
    # def on_bond_double_spin_box_editingFinished(self):
    #     self.model.bond = self.bond_double_spin_box.value()

    @pyqtSlot(float)
    def on_bond_double_spin_box_valueChanged(self, value):
        self.model.bond = value


class SWNTViewMixin:
    """Mixin class for nanotube classes."""
    # @pyqtSlot()
    # def on_swnt_n_spin_box_editingFinished(self):
    #     self.model.n = self.swnt_n_spin_box.value()

    @pyqtSlot(int)
    def on_swnt_n_spin_box_valueChanged(self, value):
        self.model.n = value

    # @pyqtSlot()
    # def on_swnt_m_spin_box_editingFinished(self):
    #     self.model.m = self.swnt_m_spin_box.value()

    @pyqtSlot(int)
    def on_swnt_m_spin_box_valueChanged(self, value):
        self.model.m = value

    # @pyqtSlot()
    # def on_swnt_nz_double_spin_box_editingFinished(self):
    #     self.model.nz = self.swnt_nz_double_spin_box.value()
    #     self._update_Lz_double_spin_box()

    @pyqtSlot(float)
    def on_swnt_nz_double_spin_box_valueChanged(self, value):
        print('on_swnt_nz_double_spin_box_valueChanged')
        self.model.nz = value

    # @pyqtSlot()
    # def on_swnt_Lz_double_spin_box_editingFinished(self):
    #     self.model.Lz = self.swnt_Lz_double_spin_box.value()

    @pyqtSlot(float)
    def on_swnt_Lz_double_spin_box_valueChanged(self, value):
        print('on_swnt_Lz_double_spin_box_valueChanged')
        self.model.Lz = value

    @pyqtSlot(int)
    def on_swnt_fix_Lz_check_box_stateChanged(self, value):
        if value:
            self.model.fix_Lz = True
            self.swnt_nz_double_spin_box.setReadOnly(True)
            self.swnt_Lz_double_spin_box.setReadOnly(False)
        else:
            self.model.fix_Lz = False
            self.swnt_nz_double_spin_box.setReadOnly(False)
            self.swnt_Lz_double_spin_box.setReadOnly(True)


class MWNTViewMixin:

    # @pyqtSlot()
    # def on_mwnt_Lz_spin_box_editingFinished(self):
    #     self.model.Lz = self.mwnt_Lz_double_spin_box.value()

    @pyqtSlot(float)
    def on_mwnt_Lz_double_spin_box_valueChanged(self, value):
        self.model.Lz = value

    @pyqtSlot(int)
    def on_mwnt_generator_button_group_buttonClicked(self, value):
        print(value)


class BundleViewMixin:

    @pyqtSlot(int)
    def on_bundle_generator_check_box_stateChanged(self, value):
        [spin_box.setReadOnly(False if value else True)
         for spin_box in (self.bundle_nx_spin_box, self.bundle_ny_spin_box)]

    # @pyqtSlot()
    # def on_bundle_nx_spin_box_editingFinished(self):
    #     self.model.nx = self.bundle_nx_spin_box.value()

    @pyqtSlot(int)
    def on_bundle_nx_spin_box_valueChanged(self, value):
        self.model.nx = value

    # @pyqtSlot()
    # def on_bundle_ny_spin_box_editingFinished(self):
    #     self.model.ny = self.bundle_ny_spin_box.value()

    @pyqtSlot(int)
    def on_bundle_ny_spin_box_valueChanged(self, value):
        self.model.ny = value


class GrapheneViewMixin:
    """Mixin class for graphene."""

    @pyqtSlot(int)
    def on_graphene_generator_button_group_buttonClicked(self, value):
        use_edge_lengths = self.achiral_edge_lengths_radio_button.isChecked()
        [obj.setReadOnly(use_edge_lengths) for obj in
         (self.swnt_n_spin_box, self.swnt_m_spin_box,
          self.unrolled_swnt_nx_spin_box,
          self.unrolled_swnt_nz_spin_box,
          self.unrolled_swnt_Lx_double_spin_box,
          self.unrolled_swnt_Lz_double_spin_box)]
        [obj.setReadOnly(not use_edge_lengths) for obj in
         (self.armchair_edge_length_double_spin_box,
          self.zigzag_edge_length_double_spin_box)]

    @pyqtSlot(int)
    def on_stacking_order_button_group_buttonClicked(self, value):
        print(value)
        # if self.nlayer_spin_box.value() == 1:
        #     pass

    # @pyqtSlot()
    # def on_armchair_edge_length_double_spin_box_editingFinished(self):
    #     self.model.armchair_edge_length = \
    #         self.armchair_edge_length_double_spin_box.value()

    @pyqtSlot(float)
    def on_armchair_edge_length_double_spin_box_valueChanged(self, value):
        self.model.armchair_edge_length = value

    # @pyqtSlot()
    # def on_zigzag_edge_length_double_spin_box_editingFinished(self):
    #     self.model.zigzag_edge_length = \
    #         self.zigzag_edge_length_double_spin_box.value()

    @pyqtSlot(float)
    def on_zigzag_edge_length_double_spin_box_valueChanged(self, value):
        self.model.zigzag_edge_length = value

    # @pyqtSlot()
    # def on_unrolled_swnt_nx_spin_box_editingFinished(self):
    #     self.model.nx = self.unrolled_swnt_nx_spin_box.value()

    @pyqtSlot(int)
    def on_unrolled_swnt_nx_spin_box_valueChanged(self, value):
        self.model.nx = value

    # @pyqtSlot()
    # def on_unrolled_swnt_Lx_double_spin_box_editingFinished(self):
    #     self.model.Lx = self.unrolled_swnt_Lx_double_spin_box.value()

    @pyqtSlot(float)
    def on_unrolled_swnt_Lx_double_spin_box_valueChanged(self, value):
        self.model.Lx = value

    # @pyqtSlot()
    # def on_unrolled_swnt_nz_spin_box_editingFinished(self):
    #     self.model.nz = self.unrolled_swnt_nz_spin_box.value()

    @pyqtSlot(int)
    def on_unrolled_swnt_nz_spin_box_valueChanged(self, value):
        self.model.nz = value

    # @pyqtSlot()
    # def on_unrolled_swnt_Lz_double_spin_box_editingFinished(self):
    #     self.model.Lz = self.unrolled_swnt_Lz_double_spin_box.value()

    @pyqtSlot(float)
    def on_unrolled_swnt_Lz_double_spin_box_valueChanged(self, value):
        self.model.Lz = value

    # @pyqtSlot()
    # def on_nlayers_spin_box_editingFinished(self):
    #     self.model.nlayers = self.nlayers_spin_box.value()
    #     #self.update_stacking_order_buttonGroup()
    #     self._set_layer_rotation_increment_double_spin_box_read_only_state()

    @pyqtSlot(int)
    def on_nlayers_spin_box_valueChanged(self, value):
        self.model.nlayers = value
        self._update_layer_rotation_increment_read_only_state()

    def _update_layer_rotation_increment_read_only_state(self):
        if self.model.nlayers == 1:
            self.layer_rotation_increment_double_spin_box.setValue(0.0)
        self.layer_rotation_increment_double_spin_box.setReadOnly(
            False if self.model.nlayers > 1 else True)
        # self.update_stacking_order_buttonGroup()

    # def update_stacking_order_buttonGroup(self):
    #     if self.nlayers_spin_box.value() == 1:
    #         for rb in (self.AA_stacking_radio_button,
    #                    self.AB_stacking_radio_button):
    #             rb.setCheckable(False)
    #     else:
    #         for rb in (self.AA_stacking_radio_button,
    #                    self.AB_stacking_radio_button):
    #             rb.setCheckable(True)

    # @pyqtSlot()
    # def on_layer_rotation_increment_double_spin_box_editingFinished(self):
    #     self.model.layer_rotation_increment = \
    #         self.layer_rotation_increment_double_spin_box.value()

    @pyqtSlot(float)
    def on_layer_rotation_increment_double_spin_box_valueChanged(self, value):
        self.model.layer_rotation_increment = value


class FullereneViewMixin:
    pass


class BulkStructureViewMixin:
    pass
