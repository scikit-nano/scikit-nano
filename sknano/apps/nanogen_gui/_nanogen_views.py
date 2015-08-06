# -*- coding: utf-8 -*-
"""
=======================================================================
NanoGen model views (:mod:`sknano.apps.nanogen_gui._nanogen_views`)
=======================================================================

.. currentmodule:: sknano.apps.nanogen_gui._nanogen_views

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import importlib
import os
import time

try:
    import PyQt5 as PyQt
except ImportError:
    try:
        import PyQt4 as PyQt
    except ImportError:
        PyQt = None

if PyQt is not None:
    try:
        from PyQt5.QtCore import pyqtSlot
        from PyQt5.QtWidgets import QMainWindow
        from ._pyqt5_ui_nanogen import Ui_NanoGen
        from ._pyqt5_ui_swnt_generator import Ui_SWNTGenerator
        from ._pyqt5_ui_mwnt_generator import Ui_MWNTGenerator
        from ._pyqt5_ui_graphene_generator import Ui_GrapheneGenerator
        from ._pyqt5_ui_fullerene_generator import Ui_FullereneGenerator
        from ._pyqt5_ui_bulk_structure_generator import \
            Ui_BulkStructureGenerator
    except ImportError:
        from PyQt4.QtCore import pyqtSlot
        from PyQt4.QtGui import QMainWindow
        from ._pyqt4_ui_nanogen import Ui_NanoGen
        from ._pyqt4_ui_swnt_generator import Ui_SWNTGenerator
        from ._pyqt4_ui_mwnt_generator import Ui_MWNTGenerator
        from ._pyqt4_ui_graphene_generator import Ui_GrapheneGenerator
        from ._pyqt4_ui_fullerene_generator import Ui_FullereneGenerator
        from ._pyqt4_ui_bulk_structure_generator import \
            Ui_BulkStructureGenerator

from sknano.core import get_fpath
from sknano.structures import get_chiral_indices_from_str
# from ._nanogen_controllers import NanoGenController
# from ._nanogen_models import NanoGenModel
from ._view_mixins import NanoGenViewMixin, SWNTViewMixin, MWNTViewMixin, \
    BundleViewMixin, GrapheneViewMixin, FullereneViewMixin, \
    BulkStructureViewMixin

__all__ = ['NanoGenView']


class NanoGenView(QMainWindow, Ui_NanoGen):
    """:mod:`~sknano.apps.nanogen_gui` MVC view class.

    Parameters
    ----------
    controller : :class:`NanoGenController`
        An instance of :class:`NanoGenController`.
    model : :class:`NanoGenModel`
        An instance of :class:`NanoGenModel`.

    """
    def __init__(self, controller=None, model=None):
        self.controller = controller
        self.model = model
        model.register_observer(self)
        super().__init__()
        self.setupUi(self)
        self.generator_controller = None

    @pyqtSlot()
    def on_swnt_generator_push_button_clicked(self):
        from ._nanogen_controllers import SWNTGeneratorController
        from ._nanogen_models import SWNTModel
        self.generator_controller = \
            SWNTGeneratorController(model=SWNTModel(), parent=self)
        self.update_app_view()

    @pyqtSlot()
    def on_mwnt_generator_push_button_clicked(self):
        from ._nanogen_controllers import MWNTGeneratorController
        from ._nanogen_models import MWNTModel
        self.generator_controller = \
            MWNTGeneratorController(model=MWNTModel(), parent=self)
        self.update_app_view()

    @pyqtSlot()
    def on_graphene_generator_push_button_clicked(self):
        from ._nanogen_controllers import GrapheneGeneratorController
        from ._nanogen_models import GrapheneModel
        self.generator_controller = \
            GrapheneGeneratorController(model=GrapheneModel(), parent=self)
        self.update_app_view()

    @pyqtSlot()
    def on_fullerene_generator_push_button_clicked(self):
        from ._nanogen_controllers import FullereneGeneratorController
        from ._nanogen_models import FullereneModel
        self.generator_controller = \
            FullereneGeneratorController(model=FullereneModel(), parent=self)
        self.update_app_view()

    @pyqtSlot()
    def on_bulk_structure_generator_push_button_clicked(self):
        from ._nanogen_controllers import BulkStructureGeneratorController
        from ._nanogen_models import BulkStructureModel
        self.generator_controller = \
            BulkStructureGeneratorController(model=BulkStructureModel(),
                                             parent=self)
        self.update_app_view()

    @pyqtSlot()
    def on_generate_push_button_clicked(self):
        if self.generator_controller is not None:

            kwargs = self.generator_controller.get_generator_parameters()
            generator_class = kwargs['generator_class']
            del kwargs['generator_class']

            structure_format = \
                str(self.structure_format_combo_box.itemText(
                    self.structure_format_combo_box.currentIndex()))

            generator = getattr(importlib.import_module('sknano.generators'),
                                generator_class)
            self.nanogen_status_bar.showMessage('Generating structure...')
            outpath, fname = os.path.split(self.fpath_line_edit.text())
            generator(**kwargs).save(fname=fname, outpath=outpath,
                                     structure_format=structure_format)
            self.nanogen_status_bar.showMessage('Done!')
            time.sleep(2)
            self.nanogen_status_bar.showMessage('Ready.')
        self.update_app_view()

    # @pyqtSlot()
    # def on_save_as_push_button_clicked(self):
    #     dialog = QFileDialog()

    def update_app_view(self):
        self.nanogen_status_bar.showMessage('Ready.')
        if self.generator_controller is not None:
            kwargs = self.generator_controller.get_generator_parameters()
            generator = getattr(importlib.import_module('sknano.generators'),
                                kwargs['generator_class'])
            structure_format = \
                str(self.structure_format_combo_box.itemText(
                    self.structure_format_combo_box.currentIndex()))

            fpath = get_fpath(fname=generator.generate_fname(**kwargs),
                              ext=structure_format, outpath=os.getcwd(),
                              overwrite=False, add_fnum=True)
            self.fpath_line_edit.setText(fpath)

        # self.Ch_line_edit.setText('{:.3f}'.format(self.model.Ch))
        # self.dt_line_edit.setText('{:.3f}'.format(self.model.dt))
        # self.T_line_edit.setText('{:.3f}'.format(self.model.T))
        # self.chiral_angle_line_edit.setText(
        #    '{:.2f}'.format(self.model.chiral_angle))
        # self.N_line_edit.setText(str(self.model.N))
        # self.Natoms_line_edit.setText(str(self.model.Natoms))
        # self.tube_mass_line_edit.setText(
        #    '{:.3e}'.format(self.model.tube_mass))
        # self.Natoms_per_tube_line_edit.setText(
        #    str(self.model.Natoms_per_tube))
        # self.Ntubes_mantissa_spin_box.setValue(
        #    str(self.model.Ntubes))

        # self.d_line_edit.setText(str(self.model.d))
        # self.dR_line_edit.setText(str(self.model.dR))
        # self.t1_line_edit.setText(str(self.model.t1))
        # self.t2_line_edit.setText(str(self.model.t2))


class SWNTGeneratorView(QMainWindow, Ui_SWNTGenerator, NanoGenViewMixin,
                        SWNTViewMixin, BundleViewMixin):

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

        kwargs['n'] = self.swnt_n_spin_box.value()
        kwargs['m'] = self.swnt_m_spin_box.value()
        kwargs['fix_Lz'] = self.swnt_fix_Lz_check_box.isChecked()
        kwargs['Lz'] = self.swnt_Lz_double_spin_box.value()
        kwargs['nz'] = self.swnt_nz_double_spin_box.value()

        kwargs['generator_class'] = 'SWNTGenerator'
        if self.bundle_generator_check_box.isChecked():
            kwargs['bundle_packing'] = 'hcp'
            kwargs['nx'] = self.bundle_nx_spin_box.value()
            kwargs['ny'] = self.bundle_ny_spin_box.value()
            # kwargs['Ntubes'] = self.model.structure.Ntubes
            kwargs['generator_class'] = 'SWNTBundleGenerator'

        return kwargs

    def update_app_view(self):
        self.elements_bond_label.setText('-'.join((self.model.element1,
                                                   self.model.element2,
                                                   ' bond =')))
        self.bond_double_spin_box.setValue(self.model.bond)

        self.swnt_n_spin_box.setValue(self.model.n)
        self.swnt_m_spin_box.setValue(self.model.m)

        self.swnt_nz_double_spin_box.blockSignals(True)
        self.swnt_Lz_double_spin_box.blockSignals(True)
        self.swnt_Lz_double_spin_box.setValue(self.model.Lz)
        self.swnt_nz_double_spin_box.setValue(self.model.nz)
        self.swnt_nz_double_spin_box.blockSignals(False)
        self.swnt_Lz_double_spin_box.blockSignals(False)

        self.bundle_nx_spin_box.setValue(self.model.nx)
        self.bundle_ny_spin_box.setValue(self.model.ny)

        self.bundle_Lx_line_edit.setText('{:.4f} nm'.format(self.model.Lx))
        self.bundle_Ly_line_edit.setText('{:.4f} nm'.format(self.model.Ly))
        self.parent.update_app_view()


class MWNTGeneratorView(QMainWindow, Ui_MWNTGenerator, NanoGenViewMixin,
                        MWNTViewMixin, BundleViewMixin):
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
        kwargs['Lz'] = self.mwnt_Lz_double_spin_box.value()
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
            kwargs['nx'] = self.bundle_nx_spin_box.value()
            kwargs['ny'] = self.bundle_ny_spin_box.value()
            # kwargs['Ntubes'] = self.model.structure.Ntubes
            kwargs['generator_class'] = 'MWNTBundleGenerator'
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
        self.mwnt_Lz_double_spin_box.setValue(self.model.Lz)

        self.bundle_nx_spin_box.setValue(self.model.nx)
        self.bundle_ny_spin_box.setValue(self.model.ny)
        try:
            self.bundle_Lx_line_edit.setText('{:.4f} nm'.format(self.model.Lx))
            self.bundle_Ly_line_edit.setText('{:.4f} nm'.format(self.model.Ly))
        except IndexError:
            self.bundle_Lx_line_edit.setText('{:.4f} nm'.format(0.0))
            self.bundle_Ly_line_edit.setText('{:.4f} nm'.format(0.0))
        self.parent.update_app_view()


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
            kwargs['nx'] = self.unrolled_swnt_nx_spin_box.value()
            kwargs['nz'] = self.unrolled_swnt_nz_spin_box.value()
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


class FullereneGeneratorView(QMainWindow, Ui_FullereneGenerator,
                             FullereneViewMixin):
    def __init__(self, controller=None, model=None, parent=None):
        self.controller = controller
        self.model = model
        self.parent = parent
        model.register_observer(self)
        super().__init__(parent=parent)
        self.setupUi(self)

    def get_generator_parameters(self):
        kwargs = {}
        return kwargs

    def update_app_view(self):
        self.fullerene_list_widget.addItems(
            [fullerene for fullerene in self.model.fullerenes])
        self.parent.update_app_view()


class BulkStructureGeneratorView(QMainWindow, Ui_BulkStructureGenerator,
                                 BulkStructureViewMixin):
    def __init__(self, controller=None, model=None, parent=None):
        self.controller = controller
        self.model = model
        self.parent = parent
        model.register_observer(self)
        super().__init__(parent=parent)
        self.setupUi(self)

    def get_generator_parameters(self):
        kwargs = {}
        return kwargs

    def update_app_view(self):
        self.parent.update_app_view()
