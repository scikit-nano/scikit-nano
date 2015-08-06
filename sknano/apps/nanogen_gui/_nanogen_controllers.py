# -*- coding: utf-8 -*-
"""
===============================================================================
NanoGen view controllers (:mod:`sknano.apps.nanogen_gui._nanogen_controllers`)
===============================================================================

.. currentmodule:: sknano.apps.nanogen_gui._nanogen_controllers

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'


try:
    from PyQt5.QtWidgets import QApplication
except ImportError:
    from PyQt4.QtGui import QApplication

from ._nanogen_views import NanoGenView, SWNTGeneratorView, \
    MWNTGeneratorView, GrapheneGeneratorView, \
    FullereneGeneratorView, BulkStructureGeneratorView

__all__ = ['NanoGenController', 'SWNTGeneratorController',
           'MWNTGeneratorController', 'GrapheneGeneratorController',
           'BulkStructureGeneratorController', 'ViewControllerMixin',
           'GeneratorViewController']


class ViewControllerMixin:
    def refresh_view(self):
        """Refresh `NanoGenView`."""
        self.view.update_app_view()


class NanoGenController(ViewControllerMixin):
    """:mod:`~sknano.apps.nanogen_gui` MVC controller class.

    Parameters
    ----------
    args : :attr:`python:sys.argv`
    model : :class:`~sknano.apps.nanogen_gui._nanogen_model.NanoGenModel`
        An instance of
        :class:`~sknano.apps.nanogen_gui._nanogen_model.NanoGenModel`.

    """
    def __init__(self, args, model=None):
        app = QApplication(args)
        self.model = model
        self.view = NanoGenView(self, self.model)
        self.model.notify_observers()
        self.view.show()
        app.exec_()


class GeneratorViewController(ViewControllerMixin):
    def __init__(self, model=None, view=None, **kwargs):
        self.model = model
        self.view = view(self, self.model, **kwargs)
        self.model.notify_observers()
        self.view.show()

    def get_generator_parameters(self):
        return self.view.get_generator_parameters()


class SWNTGeneratorController(GeneratorViewController):
    def __init__(self, model=None, **kwargs):
        super().__init__(model=model, view=SWNTGeneratorView, **kwargs)


class MWNTGeneratorController(GeneratorViewController):
    def __init__(self, model=None, **kwargs):
        super().__init__(model=model, view=MWNTGeneratorView, **kwargs)


class GrapheneGeneratorController(GeneratorViewController):
    def __init__(self, model=None, **kwargs):
        super().__init__(model=model, view=GrapheneGeneratorView, **kwargs)


class FullereneGeneratorController(GeneratorViewController):
    def __init__(self, model=None, **kwargs):
        super().__init__(model=model, view=FullereneGeneratorView, **kwargs)


class BulkStructureGeneratorController(GeneratorViewController):
    def __init__(self, model=None, **kwargs):
        super().__init__(model=model, view=BulkStructureGeneratorView,
                         **kwargs)
