# -*- coding: utf-8 -*-
"""
=======================================================
NanoGen view (:mod:`sknano.apps.nanogen_gui._ng_view`)
=======================================================

.. currentmodule:: sknano.apps.nanogen_gui._ng_view

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import importlib
import os
import time

try:
    # from PyQt4.QtCore import pyqtSlot
    # from PyQt4.QtGui import QMainWindow
    from PyQt5.QtCore import pyqtSlot
    from PyQt5.QtWidgets import QMainWindow
    from ._ui_nanogen import Ui_NanoGen
except ImportError as e:
    print(e)

from sknano.core import get_fpath
from ._view_mixins import MainWindowViewMixin, SWNTViewMixin, MWNTViewMixin, \
    GrapheneViewMixin, FullereneViewMixin

__all__ = ['NGView']


class NGView(QMainWindow, Ui_NanoGen, MainWindowViewMixin, SWNTViewMixin,
             MWNTViewMixin, GrapheneViewMixin, FullereneViewMixin):
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
        self.fpath = None
        self.nanogen_statusBar.showMessage('Ready.')

    @pyqtSlot()
    def on_generate_pushButton_clicked(self):
        generator_tab = \
            str(self.nanogen_tabWidget.tabText(
                self.nanogen_tabWidget.currentIndex()))
        # args = []
        kwargs = {}
        element1 = str(self.element1_comboBox.itemText(
                       self.element1_comboBox.currentIndex()))
        element2 = str(self.element2_comboBox.itemText(
                       self.element2_comboBox.currentIndex()))
        kwargs['basis'] = [element1, element2]
        kwargs['bond'] = self.bond_doubleSpinBox.value()

        if generator_tab == 'SWNTs':
            # args.extend([self.n_spinBox.value(), self.m_spinBox.value()])
            kwargs['Ch'] = (self.swnt_n_spinBox.value(),
                            self.swnt_m_spinBox.value())
            kwargs['nz'] = self.swnt_nz_doubleSpinBox.value()
            kwargs['fix_Lz'] = self.swnt_fix_Lz_checkBox.isChecked()
            kwargs['Lz'] = self.swnt_Lz_doubleSpinBox.value()

            generator_class = 'SWNTGenerator'
            if self.swnt_bundle_generator_checkBox.isChecked():
                kwargs['nx'] = self.swnt_bundle_nx_spinBox.value()
                kwargs['ny'] = self.swnt_bundle_ny_spinBox.value()
                generator_class = 'SWNTBundleGenerator'

        elif generator_tab == 'MWNTs':
            generator_class = 'MWNTGenerator'
            kwargs['Lz'] = self.mwnt_Lz_doubleSpinBox.value()

            if self.mwnt_Ch_list_radioButton.isChecked():
                pass
            else:
                kwargs['Nwalls'] = self.Nwalls_spinBox.value()
                kwargs['min_wall_diameter'] = \
                    self.min_wall_diameter_doubleSpinBox.value()
                kwargs['max_wall_diameter'] = \
                    self.max_wall_diameter_doubleSpinBox.value()
                kwargs['wall_spacing'] = \
                    self.wall_spacing_doubleSpinBox.value()

            if self.mwnt_bundle_generator_checkBox.isChecked():
                kwargs['nx'] = self.mwnt_bundle_nx_spinBox.value()
                kwargs['ny'] = self.mwnt_bundle_ny_spinBox.value()
                generator_class = 'MWNTBundleGenerator'

        elif generator_tab == 'Graphene':
            if self.achiral_edge_lengths_radioButton.isChecked():
                kwargs['armchair_edge_length'] = \
                    self.armchair_edge_length_doubleSpinBox.value()
                kwargs['zigzag_edge_length'] = \
                    self.zigzag_edge_length_doubleSpinBox.value()
                generator_class = 'GrapheneGenerator'
            else:
                kwargs['Ch'] = (self.unrolled_swnt_n_spinBox.value(),
                                self.unrolled_swnt_m_spinBox.value())
                kwargs['nx'] = self.unrolled_swnt_nx_spinBox.value()
                kwargs['nz'] = self.unrolled_swnt_nz_spinBox.value()
                generator_class = 'UnrolledSWNTGenerator'

            kwargs['layer_rotation_increment'] = \
                self.layer_rotation_increment_doubleSpinBox.value()
            kwargs['nlayers'] = self.nlayers_spinBox.value()
            # edge = 'ZZ' if self.ZZ_edge_radioButton.isChecked() else 'AC'
            kwargs['stacking_order'] = \
                'AA' if self.AA_stacking_radioButton.isChecked() else 'AB'

        fname = None
        outpath = None

        structure_format = \
            str(self.structure_format_comboBox.itemText(
                self.structure_format_comboBox.currentIndex()))

        generator = getattr(importlib.import_module('sknano.generators'),
                            generator_class)
        self.nanogen_statusBar.showMessage('Generating structure...')
        generator(**kwargs).save(fname=fname, outpath=outpath,
                                 structure_format=structure_format)
        self.nanogen_statusBar.showMessage('Done!')
        time.sleep(2)
        self.nanogen_statusBar.showMessage('Ready.')

    # @pyqtSlot()
    # def on_save_as_pushButton_clicked(self):
    #     dialog = QFileDialog()

    def _update_main_window_view(self):
        self.elements_bond_label.setText('-'.join((self.model.element1,
                                                   self.model.element2,
                                                   ' bond =')))
        self.bond_doubleSpinBox.setValue(self.model.bond)

    def _update_swnt_view(self):
        self.swnt_n_spinBox.setValue(self.model.swnt.n)
        self.swnt_m_spinBox.setValue(self.model.swnt.m)

        self.swnt_nz_doubleSpinBox.setValue(self.model.swnt.nz)
        self.swnt_Lz_doubleSpinBox.setValue(self.model.swnt.Lz)

        self.swnt_bundle_nx_spinBox.setValue(self.model.swnt_bundle.nx)
        self.swnt_bundle_ny_spinBox.setValue(self.model.swnt_bundle.ny)

        self.swnt_bundle_Lx_lineEdit.setText(
            '{:.4f} nm'.format(self.model.swnt_bundle.Lx))
        self.swnt_bundle_Ly_lineEdit.setText(
            '{:.4f} nm'.format(self.model.swnt_bundle.Ly))

    def _update_mwnt_view(self):
        self.Nwalls_spinBox.setValue(self.model.mwnt.Nwalls)
        self.min_wall_diameter_doubleSpinBox.setValue(
            self.model.mwnt.min_wall_diameter)
        self.max_wall_diameter_doubleSpinBox.setValue(
            self.model.mwnt.max_wall_diameter)
        self.wall_spacing_doubleSpinBox.setValue(self.model.mwnt.wall_spacing)
        self.mwnt_Lz_doubleSpinBox.setValue(self.model.mwnt.Lz)

        self.mwnt_bundle_nx_spinBox.setValue(self.model.mwnt_bundle.nx)
        self.mwnt_bundle_ny_spinBox.setValue(self.model.mwnt_bundle.ny)
        self.mwnt_bundle_Lx_lineEdit.setText(
            '{:.4f} nm'.format(self.model.mwnt_bundle.Lx))
        self.mwnt_bundle_Ly_lineEdit.setText(
            '{:.4f} nm'.format(self.model.mwnt_bundle.Ly))

    def _update_graphene_view(self):
        self.armchair_edge_length_doubleSpinBox.setValue(
            self.model.graphene.armchair_edge_length)
        self.zigzag_edge_length_doubleSpinBox.setValue(
            self.model.graphene.zigzag_edge_length)
        self.nlayers_spinBox.setValue(self.model.nlayers)
        self.layer_rotation_increment_doubleSpinBox(
            self.model.layer_rotation_increment)

        self.unrolled_swnt_n_spinBox.setValue(self.model.unrolled_swnt.n)
        self.unrolled_swnt_m_spinBox.setValue(self.model.unrolled_swnt.m)

    def _update_fullerene_view(self):
        pass

    def update_app_view(self):
        self._update_main_window_view()

        generator_tab = \
            str(self.nanogen_tabWidget.tabText(
                self.nanogen_tabWidget.currentIndex()))

        if generator_tab == 'SWNTs':
            self._update_swnt_view()
        elif generator_tab == 'MWNTs':
            self._update_mwnt_view()
        elif generator_tab == 'Graphene':
            self._update_graphene_view()
        elif generator_tab == 'Fullerenes':
            self._update_fullerene_view()

        if self.fpath is None:
            # args = []
            kwargs = {}

            element1 = \
                str(self.element1_comboBox.itemText(
                    self.element1_comboBox.currentIndex()))
            element2 = \
                str(self.element2_comboBox.itemText(
                    self.element2_comboBox.currentIndex()))
            kwargs['basis'] = [element1, element2]

            if generator_tab == 'SWNTs':
                kwargs['n'] = self.swnt_n_spinBox.value()
                kwargs['m'] = self.swnt_m_spinBox.value()
                if self.swnt_bundle_generator_checkBox.isChecked():
                    kwargs['nx'] = self.swnt_bundle_nx_spinBox.value()
                    kwargs['ny'] = self.swnt_bundle_ny_spinBox.value()
                    generator_class = 'SWNTBundleGenerator'
                else:
                    generator_class = 'SWNTGenerator'
                kwargs['nz'] = self.swnt_nz_doubleSpinBox.value()
            elif generator_tab == 'MWNTs':
                if self.mwnt_bundle_generator_checkBox.isChecked():
                    kwargs['nx'] = self.mwnt_bundle_nx_spinBox.value()
                    kwargs['ny'] = self.mwnt_bundle_ny_spinBox.value()
                    generator_class = 'MWNTBundleGenerator'
                else:
                    generator_class = 'MWNTGenerator'
            elif generator_tab == 'Graphene':
                if self.achiral_edge_lengths_radioButton.isChecked():
                    kwargs['armchair_edge_length'] = \
                        self.armchair_edge_length_doubleSpinBox.value()
                    kwargs['zigzag_edge_length'] = \
                        self.zigzag_edge_length_doubleSpinBox.value()
                    generator_class = 'GrapheneGenerator'
                else:
                    kwargs['nx'] = self.unrolled_swnt_nx_spinBox.value()
                    kwargs['nz'] = self.unrolled_swnt_nz_spinBox.value()
                    generator_class = 'UnrolledSWNTGenerator'

                kwargs['nlayers'] = self.nlayers_spinBox.value()

            generator = getattr(importlib.import_module('sknano.generators'),
                                generator_class)
            structure_format = \
                str(self.structure_format_comboBox.itemText(
                    self.structure_format_comboBox.currentIndex()))

            self.fpath = \
                get_fpath(fname=generator.generate_fname(**kwargs),
                          ext=structure_format, outpath=os.getcwd(),
                          overwrite=False, add_fnum=True)
            self.fpath_lineEdit.setText(self.fpath)

        # self.Ch_lineEdit.setText('{:.3f}'.format(self.model.Ch))
        # self.dt_lineEdit.setText('{:.3f}'.format(self.model.dt))
        # self.T_lineEdit.setText('{:.3f}'.format(self.model.T))
        # self.chiral_angle_lineEdit.setText(
        #    '{:.2f}'.format(self.model.chiral_angle))
        # self.N_lineEdit.setText(str(self.model.N))
        # self.Natoms_lineEdit.setText(str(self.model.Natoms))
        # self.tube_mass_lineEdit.setText(
        #    '{:.3e}'.format(self.model.tube_mass))
        # self.Natoms_per_tube_lineEdit.setText(
        #    str(self.model.Natoms_per_tube))
        # self.Ntubes_mantissa_spinBox.setValue(
        #    str(self.model.Ntubes))

        # self.d_lineEdit.setText(str(self.model.d))
        # self.dR_lineEdit.setText(str(self.model.dR))
        # self.t1_lineEdit.setText(str(self.model.t1))
        # self.t2_lineEdit.setText(str(self.model.t2))
