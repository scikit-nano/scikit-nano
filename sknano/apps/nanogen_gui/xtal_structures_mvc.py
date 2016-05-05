# -*- coding: utf-8 -*-
"""
===================================================================================
Crystal Structure MVC classes (:mod:`sknano.apps.nanogen_gui.xtal_structures_mvc`)
===================================================================================

.. currentmodule:: sknano.apps.nanogen_gui.xtal_structures_mvc

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

try:
    # from PyQt5.QtCore import pyqtSlot
    from PyQt5.QtWidgets import QMainWindow
    from ._pyqt5_ui_xtal_structure_generator import \
        Ui_CrystalStructureGenerator
except ImportError:
    # from PyQt4.QtCore import pyqtSlot
    from PyQt4.QtGui import QMainWindow
    from ._pyqt4_ui_xtal_structure_generator import \
        Ui_CrystalStructureGenerator

from .mvc_base import GeneratorViewController
from .mvc_mixins import ObserverModelMixin

__all__ = ['CrystalStructureModel', 'CrystalStructureGeneratorView',
           'CrystalStructureGeneratorController']


class CrystalStructureModel(ObserverModelMixin):
    pass


class CrystalStructureGeneratorView(QMainWindow, Ui_CrystalStructureGenerator):
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


class CrystalStructureGeneratorController(GeneratorViewController):
    def __init__(self, model=None, **kwargs):
        super().__init__(model=model, view=CrystalStructureGeneratorView,
                         **kwargs)
