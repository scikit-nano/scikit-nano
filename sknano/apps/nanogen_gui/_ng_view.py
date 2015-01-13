# -*- coding: utf-8 -*-
"""
=======================================================
NanoGen view (:mod:`sknano.apps.nanogen_gui._ng_view`)
=======================================================

.. currentmodule:: sknano.apps.nanogen_gui._ng_view

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import importlib

try:
    from PyQt4.QtCore import pyqtSlot
    from PyQt4.QtGui import QMainWindow
    from ._ui_nanogen import Ui_NanoGen
except ImportError as e:
    print(e)

from ._view_mixins import NanotubeViewMixin, GrapheneViewMixin

__all__ = ['NGView']


class NGView(QMainWindow, Ui_NanoGen, NanotubeViewMixin, GrapheneViewMixin):
    """:mod:`~sknano.apps.nanogen_gui` MVC view class.

    Parameters
    ----------
    controller : :class:`~sknano.apps.nanogen_gui._ng_controller.NGController`
        An instance of
        :class:`~sknano.apps.nanogen_gui._ng_controller.NGController`.
    model : :class:`~sknano.apps.nanogen_gui._ng_model.NGModel`
        An instance of :class:`~sknano.apps.nanogen_gui._ng_model.NGModel`.

    """
    def __init__(self, controller=None, model=None):
        self.controller = controller
        self.model = model
        model.register_observer(self)
        super().__init__()
        self.setupUi(self)
        self.nanogen_statusBar.showMessage('Ready.')

    @pyqtSlot()
    def on_generate_pushButton_clicked(self):
        generator_tab = \
            str(self.nanogen_tabWidget.tabText(
                self.nanogen_tabWidget.currentIndex()))
        kwargs = {}

        if generator_tab == 'Nanotubes':
            kwargs['n'] = self.n_spinBox.value()
            kwargs['m'] = self.m_spinBox.value()
            kwargs['nz'] = self.nz_spinBox.value()
            kwargs['bond'] = self.nanotube_bond_doubleSpinBox.value()
            if self.nanotube_generator_radioButton.isChecked():
                if self.bundle_generator_checkBox.isChecked():
                    kwargs['nx'] = self.nx_spinBox.value()
                    kwargs['ny'] = self.ny_spinBox.value()
                    generator_class = 'SWNTBundleGenerator'
                else:
                    generator_class = 'SWNTGenerator'
            elif self.mwnt_generator_radioButton.isChecked():
                if self.bundle_generator_checkBox.isChecked():
                    kwargs['nx'] = self.nx_spinBox.value()
                    kwargs['ny'] = self.ny_spinBox.value()
                    generator_class = 'MWNTBundleGenerator'
                else:
                    generator_class = 'MWNTGenerator'
            elif self.unrolled_nanotube_generator_radioButton.isChecked():
                generator_class = 'UnrolledNanotubeGenerator'
        elif generator_tab == 'Graphene':
            kwargs['width'] = self.width_doubleSpinBox.value()
            kwargs['length'] = self.length_doubleSpinBox.value()
            kwargs['bond'] = self.graphene_bond_doubleSpinBox.value()
            kwargs['nlayers'] = self.nlayers_spinBox.value()

            edge = 'AC'
            if self.ZZ_edge_radioButton.isChecked():
                edge = 'ZZ'
            kwargs['edge'] = edge

            stacking_order = 'AB'
            if self.AA_stacking_radioButton.isChecked():
                stacking_order = 'AA'
            kwargs['stacking_order'] = stacking_order

            generator_class = 'GrapheneGenerator'

        structure_format = \
            str(self.structure_format_comboBox.itemText(
                self.structure_format_comboBox.currentIndex()))
        if structure_format.endswith('data'):
            structure_format = 'data'

        generator = getattr(importlib.import_module('sknano.generators'),
                            generator_class)
        #generator(**kwargs).save_data(fname=fname,
        #                              structure_format=structure_format)
        generator(**kwargs).save_data(structure_format=structure_format)

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
