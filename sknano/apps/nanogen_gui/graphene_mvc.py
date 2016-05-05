# -*- coding: utf-8 -*-
"""
===============================================================================
Graphene MVC classes (:mod:`sknano.apps.nanogen_gui.graphene_mvc`)
===============================================================================

.. currentmodule:: sknano.apps.nanogen_gui.graphene_mvc

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

try:
    from PyQt5.QtCore import pyqtSlot
    from PyQt5.QtWidgets import QMainWindow
    from ._pyqt5_ui_graphene_generator import Ui_GrapheneGenerator
except ImportError:
    from PyQt4.QtCore import pyqtSlot
    from PyQt4.QtGui import QMainWindow
    from ._pyqt4_ui_graphene_generator import Ui_GrapheneGenerator


from sknano.core.structures import Graphene, UnrolledSWNT
from .unrolled_swnt_mvc import UnrolledSWNTModelMixin

from .mvc_base import NanoStructureModelBase, GeneratorViewController

__all__ = ['GrapheneModel', 'GrapheneViewMixin', 'GrapheneGeneratorController']


class GrapheneModel(NanoStructureModelBase, UnrolledSWNTModelMixin):

    def __init__(self):
        super().__init__()

        self.conventional_cell_graphene = \
            Graphene.from_conventional_cell(armchair_edge_length=1,
                                            zigzag_edge_length=1,
                                            basis=self.basis,
                                            bond=self.bond)
        self.primitive_cell_graphene = \
            Graphene.from_primitive_cell(edge_length=1,
                                         basis=self.basis,
                                         bond=self.bond)

        self.structure = UnrolledSWNT((10, 10), basis=self.basis,
                                      bond=self.bond, n1=1, n3=1)
        self.nlayers = 1
        self.layer_rotation_increment = 0.0
        self.notify_observers()

    @property
    def armchair_edge_length(self):
        return self.conventional_cell_graphene.armchair_edge_length

    @armchair_edge_length.setter
    def armchair_edge_length(self, value):
        self.conventional_cell_graphene.armchair_edge_length = value
        self.notify_observers()

    @property
    def zigzag_edge_length(self):
        return self.conventional_cell_graphene.zigzag_edge_length

    @zigzag_edge_length.setter
    def zigzag_edge_length(self, value):
        self.conventional_cell_graphene.zigzag_edge_length = value
        self.notify_observers()

    @property
    def edge_length(self):
        return self.primitive_cell_graphene.edge_length

    @edge_length.setter
    def edge_length(self, value):
        self.primitive_cell_graphene.edge_length = value
        self.notify_observers()

    @property
    def nlayers(self):
        return self._nlayers

    @nlayers.setter
    def nlayers(self, value):
        self._nlayers = value
        self.notify_observers()

    @property
    def layer_rotation_increment(self):
        return self._layer_rotation_increment

    @layer_rotation_increment.setter
    def layer_rotation_increment(self, value):
        self._layer_rotation_increment = value
        self.notify_observers()


class GrapheneViewMixin:
    """Mixin class for graphene."""

    @pyqtSlot(int)
    def on_graphene_generator_button_group_buttonClicked(self, value):
        use_unrolled_swnt = True if \
            self.unrolled_swnt_chirality_radio_button.isChecked() else False
        use_conventional_unit_cell = True if \
            self.conventional_unit_cell_radio_button.isChecked() else False
        use_primitive_unit_cell = True if \
            self.primitive_unit_cell_radio_button.isChecked() else False
        [obj.setReadOnly(not use_unrolled_swnt) for obj in
         (self.swnt_n_spin_box, self.swnt_m_spin_box,
          self.unrolled_swnt_n1_spin_box, self.unrolled_swnt_n3_spin_box)]
        [obj.setReadOnly(not use_conventional_unit_cell) for
         obj in (self.armchair_edge_length_double_spin_box,
                 self.zigzag_edge_length_double_spin_box)]
        self.edge_length_double_spin_box.setReadOnly(
            not use_primitive_unit_cell)
        self.update_app_view()

    @pyqtSlot(int)
    def on_stacking_order_button_group_buttonClicked(self, value):
        print(value)
        # if self.nlayer_spin_box.value() == 1:
        #     pass

    @pyqtSlot(float)
    def on_armchair_edge_length_double_spin_box_valueChanged(self, value):
        self.model.armchair_edge_length = value

    @pyqtSlot(float)
    def on_zigzag_edge_length_double_spin_box_valueChanged(self, value):
        self.model.zigzag_edge_length = value

    @pyqtSlot(float)
    def on_edge_length_double_spin_box_valueChanged(self, value):
        self.model.edge_length = value

    @pyqtSlot(int)
    def on_unrolled_swnt_n1_spin_box_valueChanged(self, value):
        self.model.n1 = value

    @pyqtSlot(int)
    def on_unrolled_swnt_n3_spin_box_valueChanged(self, value):
        self.model.n3 = value

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
        self.update_app_view()
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


class GrapheneGeneratorView(QMainWindow, Ui_GrapheneGenerator,
                            NanoGenViewMixin, GrapheneViewMixin):
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
        if self.conventional_unit_cell_radio_button.isChecked():
            kwargs['armchair_edge_length'] = \
                self.armchair_edge_length_double_spin_box.value()
            kwargs['zigzag_edge_length'] = \
                self.zigzag_edge_length_double_spin_box.value()
            kwargs['generator_class'] = 'ConventionalCellGrapheneGenerator'
        elif self.primitive_unit_cell_radio_button.isChecked():
            kwargs['edge_length'] = \
                self.edge_length_double_spin_box.value()
            kwargs['generator_class'] = 'PrimitiveCellGrapheneGenerator'
        else:
            kwargs['n'] = self.swnt_n_spin_box.value()
            kwargs['m'] = self.swnt_m_spin_box.value()
            kwargs['n1'] = self.unrolled_swnt_n1_spin_box.value()
            kwargs['n3'] = self.unrolled_swnt_n3_spin_box.value()
            kwargs['generator_class'] = 'UnrolledSWNTGenerator'

        kwargs['layer_rotation_increment'] = \
            self.layer_rotation_increment_double_spin_box.value()
        kwargs['nlayers'] = self.nlayers_spin_box.value()
        # edge = 'ZZ' if self.ZZ_edge_radio_button.isChecked() else 'AC'
        kwargs['stacking_order'] = \
            'AA' if self.AA_stacking_radio_button.isChecked() else 'AB'

        return kwargs

    def update_app_view(self):
        self.elements_bond_label.setText('-'.join((self.model.element1,
                                                   self.model.element2,
                                                   ' bond =')))
        self.bond_double_spin_box.setValue(self.model.bond)

        self.armchair_edge_length_double_spin_box.setValue(
            self.model.armchair_edge_length)
        self.zigzag_edge_length_double_spin_box.setValue(
            self.model.zigzag_edge_length)
        self.edge_length_double_spin_box.setValue(
            self.model.edge_length)
        self.nlayers_spin_box.setValue(self.model.nlayers)
        self.layer_rotation_increment_double_spin_box.setValue(
            self.model.layer_rotation_increment)

        self.swnt_n_spin_box.setValue(self.model.n)
        self.swnt_m_spin_box.setValue(self.model.m)
        self.parent.update_app_view()


class GrapheneGeneratorController(GeneratorViewController):
    def __init__(self, model=None, **kwargs):
        super().__init__(model=model, view=GrapheneGeneratorView, **kwargs)


