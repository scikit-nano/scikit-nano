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
    UnrolledNanotubeGenerator, GrapheneGenerator, BiLayerGrapheneGenerator, \
    MWNTGenerator, MWNTBundleGenerator

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

    @pyqtSlot(int)
    def on_n_spinBox_valueChanged(self, value):
        self.model.n = value

    @pyqtSlot(int)
    def on_m_spinBox_valueChanged(self, value):
        self.model.m = value

    @pyqtSlot()
    def on_bond_doubleSpinBox_editingFinished(self):
        self.model.bond = self.bond_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_bond_doubleSpinBox_valueChanged(self, value):
        self.model.bond = value

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

    @pyqtSlot()
    def on_Lz_doubleSpinBox_editingFinished(self):
        self.model.Lz = self.Lz_doubleSpinBox.value()

    @pyqtSlot(float)
    def on_Lz_doubleSpinBox_valueChanged(self, value):
        self.model.Lz = value

    @pyqtSlot()
    def on_generate_pushButton_clicked(self):
        if self.nanotube_generator_radioButton.isChecked():
            if self.bundle_generator_checkBox.isChecked():
                generator = NanotubeBundleGenerator(
                    n=self.n_spinBox.value(),
                    m=self.m_spinBox.value(),
                    nx=self.nx_spinBox.value(),
                    ny=self.ny_spinBox.value(),
                    nz=self.nz_spinBox.value(),
                    bond=self.bond_doubleSpinBox.value())
            else:
                generator = NanotubeGenerator(
                    n=self.n_spinBox.value(),
                    m=self.m_spinBox.value(),
                    nz=self.nz_spinBox.value(),
                    bond=self.bond_doubleSpinBox.value())
        elif self.mwnt_generator_radioButton.isChecked():
            if self.bundle_generator_checkBox.isChecked():
                generator = MWNTBundleGenerator(
                    n=self.n_spinBox.value(),
                    m=self.m_spinBox.value(),
                    nx=self.nx_spinBox.value(),
                    ny=self.ny_spinBox.value(),
                    nz=self.nz_spinBox.value(),
                    bond=self.bond_doubleSpinBox.value())
            else:
                generator = MWNTGenerator(
                    n=self.n_spinBox.value(),
                    m=self.m_spinBox.value(),
                    nz=self.nz_spinBox.value(),
                    bond=self.bond_doubleSpinBox.value())
        elif self.unrolled_nanotube_generator_radioButton.isChecked():
            generator = UnrolledNanotubeGenerator(
                n=self.n_spinBox.value(),
                m=self.m_spinBox.value(),
                nx=self.nx_spinBox.value(),
                nz=self.nz_spinBox.value(),
                bond=self.bond_doubleSpinBox.value())

        structure_format = \
            str(self.structure_format_comboBox.itemText(
                self.structure_format_comboBox.currentIndex()))
        if structure_format.endswith('data'):
            structure_format = 'data'
        generator.save_data(structure_format=structure_format)

    def update_app_view(self):
        self.n_spinBox.setValue(self.model.n)
        self.m_spinBox.setValue(self.model.m)
        self.bond_doubleSpinBox.setValue(self.model.bond)
        self.nx_spinBox.setValue(self.model.nx)
        self.ny_spinBox.setValue(self.model.ny)
        self.nz_spinBox.setValue(self.model.nz)
        self.Lx_doubleSpinBox.setValue(self.model.Lx)
        self.Ly_doubleSpinBox.setValue(self.model.Ly)
        self.Lz_doubleSpinBox.setValue(self.model.Lz)

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
