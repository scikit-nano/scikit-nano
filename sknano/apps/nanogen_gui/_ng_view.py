# -*- coding: utf-8 -*-
"""
==================================================
NanoGen view (:mod:`sknano.nanogen_gui._ng_view`)
==================================================

.. currentmodule:: sknano.nanogen_gui._ng_view

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from PyQt4.QtCore import pyqtSlot
from PyQt4.QtGui import QMainWindow

from ..nanogen import NanotubeGenerator, NanotubeBundleGenerator, \
    MWNTGenerator, MWNTBundleGenerator, UnrolledNanotubeGenerator, \
    GrapheneGenerator  # , BiLayerGrapheneGenerator
from ._ui_nanogen import Ui_NanoGen

__all__ = ['NGView']


class NGView(QMainWindow, Ui_NanoGen):
    """:mod:`~sknano.nanogen_gui` MVC view class.

    Parameters
    ----------
    controller : :class:`~sknano.nanogen_gui._ng_controller.NGController`
        An instance of
        :class:`~sknano.nanogen_gui._ng_controller.NGController`.
    model : :class:`~sknano.nanogen_gui._ng_model.NGModel`
        An instance of :class:`~sknano.nanogen_gui._ng_model.NGModel`.

    """
    def __init__(self, controller=None, model=None):
        self.controller = controller
        self.model = model
        model.register_observer(self)
        super(NGView, self).__init__()
        self.setupUi(self)
        self.nanogen_statusBar.showMessage('Ready')

    # Nanotube Generator slots/signals

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

    # Graphene Generator slots/signals

    #@pyqtSlot(int)
    #def on_stacking_order_buttonGroup_buttonClicked(self):
    #    if self.nlayer_spinBox.value() == 1:
    #        pass

    @pyqtSlot()
    def on_length_doubleSpinBox_editingFinished(self):
        self.model.length = self.length_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_length_doubleSpinBox_valueChanged(self, value):
        self.model.length = value

    @pyqtSlot()
    def on_width_doubleSpinBox_editingFinished(self):
        self.model.width = self.width_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_width_doubleSpinBox_valueChanged(self, value):
        self.model.width = value

    @pyqtSlot()
    def on_nlayers_spinBox_editingFinished(self):
        self.model.nlayers = self.nlayers_spinBox.value()
        #self.update_stacking_order_buttonGroup()

    @pyqtSlot(int)
    def on_nlayers_spinBox_valueChanged(self, value):
        self.model.nlayers = value
        #self.update_stacking_order_buttonGroup()

    #def update_stacking_order_buttonGroup(self):
    #    if self.nlayers_spinBox.value() == 1:
    #        for rb in (self.AA_stacking_radioButton,
    #                   self.AB_stacking_radioButton):
    #            rb.setCheckable(False)
    #    else:
    #        for rb in (self.AA_stacking_radioButton,
    #                   self.AB_stacking_radioButton):
    #            rb.setCheckable(True)

    @pyqtSlot()
    def on_graphene_bond_doubleSpinBox_editingFinished(self):
        self.model.graphene_bond = self.graphene_bond_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_graphene_bond_doubleSpinBox_valueChanged(self, value):
        self.model.graphene_bond = value

    @pyqtSlot()
    def on_generate_pushButton_clicked(self):
        generator_tab = \
            str(self.nanogen_tabWidget.tabText(
                self.nanogen_tabWidget.currentIndex()))

        if generator_tab == 'Nanotubes':
            if self.nanotube_generator_radioButton.isChecked():
                if self.bundle_generator_checkBox.isChecked():
                    generator = NanotubeBundleGenerator(
                        n=self.n_spinBox.value(),
                        m=self.m_spinBox.value(),
                        nx=self.nx_spinBox.value(),
                        ny=self.ny_spinBox.value(),
                        nz=self.nz_spinBox.value(),
                        bond=self.nanotube_bond_doubleSpinBox.value())
                else:
                    generator = NanotubeGenerator(
                        n=self.n_spinBox.value(),
                        m=self.m_spinBox.value(),
                        nz=self.nz_spinBox.value(),
                        bond=self.nanotube_bond_doubleSpinBox.value())
            elif self.mwnt_generator_radioButton.isChecked():
                if self.bundle_generator_checkBox.isChecked():
                    generator = MWNTBundleGenerator(
                        n=self.n_spinBox.value(),
                        m=self.m_spinBox.value(),
                        nx=self.nx_spinBox.value(),
                        ny=self.ny_spinBox.value(),
                        nz=self.nz_spinBox.value(),
                        bond=self.nanotube_bond_doubleSpinBox.value())
                else:
                    generator = MWNTGenerator(
                        n=self.n_spinBox.value(),
                        m=self.m_spinBox.value(),
                        nz=self.nz_spinBox.value(),
                        bond=self.nanotube_bond_doubleSpinBox.value())
            elif self.unrolled_nanotube_generator_radioButton.isChecked():
                generator = UnrolledNanotubeGenerator(
                    n=self.n_spinBox.value(),
                    m=self.m_spinBox.value(),
                    nx=self.nx_spinBox.value(),
                    nz=self.nz_spinBox.value(),
                    bond=self.nanotube_bond_doubleSpinBox.value())
        elif generator_tab == 'Graphene':
            edge = 'AC'
            if self.ZZ_edge_radioButton.isChecked():
                edge = 'ZZ'
            stacking_order = None
            if self.AA_stacking_radioButton.isChecked():
                stacking_order = 'AA'
            elif self.AB_stacking_radioButton.isChecked():
                stacking_order = 'AB'
            generator = GrapheneGenerator(
                width=self.width_doubleSpinBox.value(),
                length=self.length_doubleSpinBox.value(),
                edge=edge,
                bond=self.graphene_bond_doubleSpinBox.value(),
                nlayers=self.nlayers_spinBox.value(),
                stacking_order=stacking_order)

        structure_format = \
            str(self.structure_format_comboBox.itemText(
                self.structure_format_comboBox.currentIndex()))
        if structure_format.endswith('data'):
            structure_format = 'data'

        generator.save_data(structure_format=structure_format)

    def update_app_view(self):
        self.n_spinBox.setValue(self.model.n)
        self.m_spinBox.setValue(self.model.m)
        self.nanotube_bond_doubleSpinBox.setValue(self.model.nanotube_bond)
        self.nx_spinBox.setValue(self.model.nx)
        self.ny_spinBox.setValue(self.model.ny)
        self.nz_spinBox.setValue(self.model.nz)
        self.Lx_lineEdit.setText('{:.4f} nm'.format(self.model.Lx))
        self.Ly_lineEdit.setText('{:.4f} nm'.format(self.model.Ly))
        self.Lz_doubleSpinBox.setValue(self.model.Lz)

        self.width_doubleSpinBox.setValue(self.model.width)
        self.length_doubleSpinBox.setValue(self.model.length)
        self.nlayers_spinBox.setValue(self.model.nlayers)
        self.graphene_bond_doubleSpinBox.setValue(self.model.graphene_bond)

        nanotube_e1e2_bond_label_text = \
            '-'.join((self.model.nanotube_element1,
                      self.model.nanotube_element2)) + ' bond ='
        self.nanotube_e1e2_bond_label.setText(nanotube_e1e2_bond_label_text)

        graphene_e1e2_bond_label_text = \
            '-'.join((self.model.graphene_element1,
                      self.model.graphene_element2)) + ' bond ='
        self.graphene_e1e2_bond_label.setText(graphene_e1e2_bond_label_text)

        #self.Ch_lineEdit.setText('{:.3f}'.format(self.model.Ch))
        #self.dt_lineEdit.setText('{:.3f}'.format(self.model.dt))
        #self.T_lineEdit.setText('{:.3f}'.format(self.model.T))
        #self.chiral_angle_lineEdit.setText(
        #    '{:.2f}'.format(self.model.chiral_angle))
        #self.N_lineEdit.setText(str(self.model.N))
        #self.Natoms_lineEdit.setText(str(self.model.Natoms))
        #self.tube_mass_lineEdit.setText(
        #    '{:.3e}'.format(self.model.tube_mass))
        #self.Natoms_per_tube_lineEdit.setText(
        #    str(self.model.Natoms_per_tube))
        #self.Ntubes_mantissa_spinBox.setValue(
        #    str(self.model.Ntubes))

        #self.d_lineEdit.setText(str(self.model.d))
        #self.dR_lineEdit.setText(str(self.model.dR))
        #self.t1_lineEdit.setText(str(self.model.t1))
        #self.t2_lineEdit.setText(str(self.model.t2))
