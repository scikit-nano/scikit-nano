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

from ._ui_nanogen import Ui_NanoGen

__all__ = ['NGView']


class NGView(QMainWindow, Ui_NanoGen):
    """:py:mod:`~sknano.nanogen_gui` MVC view class.

    Parameters
    ----------
    controller : :py:class:`~sknano.nanogen_gui._ng_controller.NGController`
        An instance of
        :py:class:`~sknano.nanogen_gui._ng_controller.NGController`.
    model : :py:class:`~sknano.nanogen_gui._ng_model.NGModel`
        An instance of :py:class:`~sknano.nanogen_gui._ng_model.NGModel`.

    """
    def __init__(self, controller=None, model=None):
        self.controller = controller
        self.model = model
        model.register_observer(self)
        super(NGView, self).__init__()
        self.setupUi(self)

    #def init(self):
    #    self.show()

    @pyqtSlot(int)
    def on_n_spinBox_valueChanged(self, value):
        self.model.n = value

    @pyqtSlot(int)
    def on_m_spinBox_valueChanged(self, value):
        self.model.m = value

    #@pyqtSlot()
    #def on_bond_doubleSpinBox_editingFinished(self):
    #    self.model.bond = self.bond_doubleSpinBox.value()

    #@pyqtSlot(float)
    #def on_bond_doubleSpinBox_valueChanged(self, value):
    #    self.model.bond = value

    @pyqtSlot()
    def on_nz_spinBox_editingFinished(self):
        self.model.nz = self.nz_spinBox.value()

    @pyqtSlot(float)
    def on_nz_spinBox_valueChanged(self, value):
        self.model.nz = int(value)

    #@pyqtSlot(int)
    #def on_fix_buttonGroup_buttonClicked(self):
    #    if self.fix_Ncells_radioButton.isChecked():
    #        self.model.fix_nz = True
    #        self.model.fix_tube_length = False
    #        self.nz_doubleSpinBox.setReadOnly(True)
    #        self.tube_length_lineEdit.setReadOnly(False)
    #    else:
    #        self.model.fix_nz = False
    #        self.model.fix_tube_length = True
    #        self.nz_doubleSpinBox.setReadOnly(False)
    #        self.tube_length_lineEdit.setReadOnly(True)

    @pyqtSlot()
    def on_Lz_lineEdit_editingFinished(self):
        self.model.Lz = self.Lz_doubleSpinBox.value()

    def update_app_view(self):
        self.n_spinBox.setValue(self.model.n)
        self.m_spinBox.setValue(self.model.m)

        #self.Ch_lineEdit.setText('{:.3f}'.format(self.model.Ch))
        #self.dt_lineEdit.setText('{:.3f}'.format(self.model.dt))
        #self.T_lineEdit.setText('{:.3f}'.format(self.model.T))
        #self.chiral_angle_lineEdit.setText(
        #    '{:.2f}'.format(self.model.chiral_angle))
        #self.N_lineEdit.setText(str(self.model.N))
        #self.Natoms_lineEdit.setText(str(self.model.Natoms))
        #self.bond_doubleSpinBox.setValue(self.model.bond)

        self.nz_spinBox.setValue(self.model.nz)
        self.Lz_doubleSpinBox.setValue(self.model.Lz)
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
