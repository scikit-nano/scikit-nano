# -*- coding: utf-8 -*-
"""
===============================================================================
NanoGen MVC classes (:mod:`sknano.apps.nanogen_gui.nanogen_mvc`)
===============================================================================

.. currentmodule:: sknano.apps.nanogen_gui.nanogen_mvc

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
        from PyQt5.QtWidgets import QApplication, QMainWindow
        from ._pyqt5_ui_nanogen import Ui_NanoGen
    except ImportError:
        from PyQt4.QtCore import pyqtSlot
        from PyQt4.QtGui import QApplication, QMainWindow
        from ._pyqt4_ui_nanogen import Ui_NanoGen


from sknano.core import get_fpath
from .mvc_mixins import ObserverModelMixin, ViewControllerMixin


__all__ = ['NanoGenModel', 'NanoGenView', 'NanoGenController']


class NanoGenModel(ObserverModelMixin):
    """:mod:`~sknano.apps.nanogen_gui` MVC model class."""
    def __init__(self):
        self._observers = []
        self.notify_observers()


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
        from .swnt_mvc import SWNTGeneratorController, SWNTModel
        self.generator_controller = \
            SWNTGeneratorController(model=SWNTModel(), parent=self)
        self.update_app_view()

    @pyqtSlot()
    def on_mwnt_generator_push_button_clicked(self):
        print('in on_mwnt_generator_push_button_clicked')
        from .mwnt_mvc import MWNTGeneratorController, MWNTModel
        self.generator_controller = \
            MWNTGeneratorController(model=MWNTModel(), parent=self)
        self.update_app_view()

    @pyqtSlot()
    def on_graphene_generator_push_button_clicked(self):
        from .graphene_mvc import GrapheneGeneratorController, GrapheneModel
        self.generator_controller = \
            GrapheneGeneratorController(model=GrapheneModel(), parent=self)
        self.update_app_view()

    @pyqtSlot()
    def on_fullerene_generator_push_button_clicked(self):
        from .fullerene_mvc import FullereneGeneratorController, FullereneModel
        self.generator_controller = \
            FullereneGeneratorController(model=FullereneModel(), parent=self)
        self.update_app_view()

    @pyqtSlot()
    def on_xtal_structure_generator_push_button_clicked(self):
        from .xtal_structures_mvc import CrystalStructureGeneratorController, \
            CrystalStructureModel
        self.generator_controller = \
            CrystalStructureGeneratorController(model=CrystalStructureModel(),
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


class NanoGenController(ViewControllerMixin):
    """:mod:`~sknano.apps.nanogen_gui` MVC controller class.

    Parameters
    ----------
    args : :attr:`python:sys.argv`
    model : :class:`~sknano.apps.nanogen_gui.nanogen_model.NanoGenModel`
        An instance of
        :class:`~sknano.apps.nanogen_gui.nanogen_model.NanoGenModel`.

    """
    def __init__(self, args, model=None):
        app = QApplication(args)
        self.model = model
        self.view = NanoGenView(self, self.model)
        self.model.notify_observers()
        self.view.show()
        app.exec_()
