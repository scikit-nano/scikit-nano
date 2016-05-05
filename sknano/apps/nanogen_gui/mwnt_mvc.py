# -*- coding: utf-8 -*-
"""
===============================================================================
MWNT MVC classes (:mod:`sknano.apps.nanogen_gui.mwnt_mvc`)
===============================================================================

.. currentmodule:: sknano.apps.nanogen_gui.mwnt_mvc

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

try:
    from PyQt5.QtCore import pyqtSlot
    from PyQt5.QtWidgets import QMainWindow
    from PyQt5.QtWidgets import QDialog
    from ._pyqt5_ui_mwnt_Ch_list_item_dialog import Ui_MWNTChListItemDialog
    from ._pyqt5_ui_mwnt_generator import Ui_MWNTGenerator
except ImportError:
    from PyQt4.QtCore import pyqtSlot
    from PyQt4.QtGui import QMainWindow
    from PyQt4.QtGui import QDialog
    from ._pyqt4_ui_mwnt_Ch_list_item_dialog import Ui_MWNTChListItemDialog
    from ._pyqt4_ui_mwnt_generator import Ui_MWNTGenerator

from sknano.core.structures import MWNT, get_chiral_indices_from_str

from .mvc_base import NanoStructureModelBase, GeneratorViewController
from .mvc_mixins import NanoStructureViewMixin
from .nanotube_bundle_mvc_mixins import NanotubeBundleModelMixin, \
    NanotubeBundleViewMixin

__all__ = ['MWNTModelMixin', 'MWNTViewMixin', 'MWNTModel', 'MWNTGeneratorView',
           'MWNTGeneratorController']


class MWNTModelMixin:
    @property
    def L(self):
        return self.structure.L

    @L.setter
    def L(self, value):
        self.structure.L = value
        self.notify_observers()


class MWNTChListItemDialog(QDialog, Ui_MWNTChListItemDialog):
    def __init__(self, item=None, parent=None):
        self.item = item
        self.parent = parent
        super().__init__(parent)
        self.setupUi(self)
        if item is not None:
            n, m = get_chiral_indices_from_str(item.text())
            self.n_spin_box.setValue(n)
            self.m_spin_box.setValue(m)

    @pyqtSlot()
    def on_ok_push_button_clicked(self):
        Ch = self.n_spin_box.value(), self.m_spin_box.value()
        if self.item is None:
            self.parent.model.Ch_list.append(Ch)
            # self.parent.mwnt_Ch_list_widget.addItem(str(Ch))
        else:
            self.item.setText(str(Ch))
        self.reject()

    @pyqtSlot()
    def on_cancel_push_button_clicked(self):
        self.reject()


class MWNTViewMixin:
    """Mixin class for MWNT generator view."""
    @pyqtSlot()
    def on_add_Ch_push_button_clicked(self):
        dialog = MWNTChListItemDialog(parent=self)
        dialog.show()
        self.update_app_view()

    @pyqtSlot()
    def on_edit_selected_Ch_push_button_clicked(self):
        selection = self.mwnt_Ch_list_widget.selectedItems()
        if len(selection) == 1:
            dialog = MWNTChListItemDialog(item=selection[0], parent=self)
            dialog.show()
            self.update_app_view()

    @pyqtSlot()
    def on_remove_selected_Ch_push_button_clicked(self):
        print(self.mwnt_Ch_list_widget.currentRow())
        self.mwnt_Ch_list_widget.takeItem(
            self.mwnt_Ch_list_widget.currentRow())
        self.update_app_view()

    @pyqtSlot()
    def on_clear_Ch_list_push_button_clicked(self):
        self.model.Ch_list.clear()
        self.update_app_view()

    @pyqtSlot(float)
    def on_mwnt_L_double_spin_box_valueChanged(self, value):
        self.model.L = value

    @pyqtSlot(float)
    def on_Nwalls_double_spin_box_valueChanged(self, value):
        self.model.Nwalls = value

    @pyqtSlot(float)
    def on_min_wall_diameter_double_spin_box_valueChanged(self, value):
        self.model.min_wall_diameter = value

    @pyqtSlot(float)
    def on_max_wall_diameter_double_spin_box_valueChanged(self, value):
        self.model.max_wall_diameter = value

    @pyqtSlot(float)
    def on_wall_spacing_double_spin_box_valueChanged(self, value):
        self.model.wall_spacing = value

    # @pyqtSlot(int)
    # def on_mwnt_generator_button_group_buttonClicked(self, value):
    #     if self.mwnt_wall_parameters_radio_button.isChecked():
    #         self.model.Nwalls = self.Nwalls_spin_box.value()
    #         self.model.min_wall_diameter = \
    #             self.min_wall_diameter_double_spin_box.value()
    #         self.model.max_wall_diameter = \
    #             self.max_wall_diameter_double_spin_box.value()
    #         self.model.wall_spacing = \
    #             self.wall_spacing_double_spin_box.value()
    #     else:
    #         self.mwnt_Ch_list_widget.clear()
    #         self.mwnt_Ch_list_widget.addItems([str(Ch) for Ch in
    #                                            self.model.Ch_list])


class MWNTModel(MWNTModelMixin, NanotubeBundleModelMixin,
                NanoStructureModelBase):
    def __init__(self):
        super().__init__()
        self.structure = \
            MWNT(Ch_list=[(10, 10), (20, 20), (30, 30)],
                 basis=self.basis, bond=self.bond, n1=1, n2=1, L=5,
                 bundle_packing='hcp')

        self.notify_observers()

    @property
    def Ch(self):
        return self.Ch_list[-1]

    @property
    def Ch_list(self):
        # return self._Ch_list
        return self.structure.Ch_list

    @Ch_list.setter
    def Ch_list(self, value):
        # self._Ch_list = value
        # self.structure.Ch_list = self.Ch_list
        self.structure.Ch_list = value
        self.notify_observers()

    @property
    def chiral_types(self):
        return self.structure.chiral_types

    @chiral_types.setter
    def chiral_types(self, value):
        self.structure.chiral_types = value
        self.notify_observers()

    @property
    def Nwalls(self):
        return self.structure.Nwalls

    @Nwalls.setter
    def Nwalls(self, value):
        self.structure.Nwalls = value
        self.notify_observers()

    @property
    def min_wall_diameter(self):
        return self.structure.min_wall_diameter

    @min_wall_diameter.setter
    def min_wall_diameter(self, value):
        self.structure.min_wall_diameter = value
        self.notify_observers()

    @property
    def max_wall_diameter(self):
        return self.structure.max_wall_diameter

    @max_wall_diameter.setter
    def max_wall_diameter(self, value):
        self.structure.max_wall_diameter = value
        self.notify_observers()

    @property
    def wall_spacing(self):
        return self.structure.wall_spacing

    @wall_spacing.setter
    def wall_spacing(self, value):
        self.structure.wall_spacing = value
        self.notify_observers()

    # def generate_Ch_list(self):
    #     self.structure = \
    #         MWNT(Nwalls=self.structure.Nwalls,
    #                    min_wall_diameter=self.structure.min_wall_diameter,
    #                    max_wall_diameter=self.structure.max_wall_diameter,
    #                    wall_spacing=self.structure.wall_spacing)


class MWNTGeneratorView(QMainWindow, Ui_MWNTGenerator, NanoStructureViewMixin,
                        MWNTViewMixin, NanotubeBundleViewMixin):
    def __init__(self, controller=None, model=None, parent=None):
        self.controller = controller
        self.model = model
        self.parent = parent
        model.register_observer(self)
        super().__init__(parent=parent)
        self.setupUi(self)

    def get_generator_parameters(self):
        kwargs = {}
        kwargs['generator_class'] = 'MWNTGenerator'
        kwargs['L'] = self.mwnt_L_double_spin_box.value()
        kwargs['Ch_list'] = \
            [get_chiral_indices_from_str(i.text()) for i in
             [self.mwnt_Ch_list_widget.item(index) for
              index in range(self.mwnt_Ch_list_widget.count())]]

        if self.mwnt_Ch_list_radio_button.isChecked():
            pass
        else:
            kwargs['Nwalls'] = self.Nwalls_spin_box.value()
            kwargs['min_wall_diameter'] = \
                self.min_wall_diameter_double_spin_box.value()
            kwargs['max_wall_diameter'] = \
                self.max_wall_diameter_double_spin_box.value()
            kwargs['wall_spacing'] = \
                self.wall_spacing_double_spin_box.value()

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

        self.mwnt_Ch_list_widget.clear()
        self.mwnt_Ch_list_widget.addItems(
            [str(Ch) for Ch in self.model.Ch_list])

        self.Nwalls_spin_box.setValue(self.model.Nwalls)
        self.min_wall_diameter_double_spin_box.setValue(
            self.model.min_wall_diameter)
        self.max_wall_diameter_double_spin_box.setValue(
            self.model.max_wall_diameter)
        self.wall_spacing_double_spin_box.setValue(
            self.model.wall_spacing)
        self.mwnt_L_double_spin_box.setValue(self.model.L)

        self.bundle_n1_spin_box.setValue(self.model.n1)
        self.bundle_n2_spin_box.setValue(self.model.n2)
        # try:
        #     self.bundle_l1_line_edit.setText('{:.4f} Å'.format(self.model.l1))
        #     self.bundle_l2_line_edit.setText('{:.4f} Å'.format(self.model.l2))
        # except IndexError:
        #     self.bundle_l1_line_edit.setText('{:.4f} Å'.format(0.0))
        #     self.bundle_l2_line_edit.setText('{:.4f} Å'.format(0.0))
        self.parent.update_app_view()


class MWNTGeneratorController(GeneratorViewController):
    def __init__(self, model=None, **kwargs):
        print('in MWNTGeneratorController')
        print('kwargs: {}'.format(kwargs))
        super().__init__(model=model, view=MWNTGeneratorView, **kwargs)
