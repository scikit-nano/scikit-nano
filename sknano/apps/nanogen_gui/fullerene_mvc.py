# -*- coding: utf-8 -*-
"""
===============================================================================
Fullerene MVC classes (:mod:`sknano.apps.nanogen_gui.fullerene_mvc`)
===============================================================================

.. currentmodule:: sknano.apps.nanogen_gui.fullerene_mvc

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

try:
    from PyQt5.QtCore import pyqtSlot
    from PyQt5.QtWidgets import QMainWindow
    from ._pyqt5_ui_fullerene_generator import Ui_FullereneGenerator
except ImportError:
    from PyQt4.QtCore import pyqtSlot
    from PyQt4.QtGui import QMainWindow
    from ._pyqt4_ui_fullerene_generator import Ui_FullereneGenerator

from sknano.core.structures import Fullerene
from .mvc_base import GeneratorViewController
from .mvc_mixins import ObserverModelMixin

__all__ = ['FullereneModel', 'FullereneGeneratorView',
           'FullereneGeneratorController']


class FullereneModel(ObserverModelMixin):
    def __init__(self):
        self._observers = []
        super().__init__()
        self.structure = Fullerene()
        self.notify_observers()

    @property
    def fullerenes(self):
        return self.structure.fullerenes


class FullereneGeneratorView(QMainWindow, Ui_FullereneGenerator):
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


class FullereneGeneratorController(GeneratorViewController):
    def __init__(self, model=None, **kwargs):
        super().__init__(model=model, view=FullereneGeneratorView, **kwargs)
