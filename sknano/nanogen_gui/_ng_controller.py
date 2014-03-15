# -*- coding: utf-8 -*-
"""
==============================================================
NanoGen controller (:mod:`sknano.nanogen_gui._ng_controller`)
==============================================================

.. currentmodule:: sknano.nanogen_gui._ng_controller

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from PyQt4.QtGui import QApplication

from ._ng_view import NGView

__all__ = ['NGController']


class NGController(object):
    """:mod:`~sknano.nanogen_gui` MVC controller class.

    Parameters
    ----------
    args : :attr:`python:sys.argv`
    model : :class:`~sknano.nanogen_gui._ng_model.NGModel`
        An instance of :class:`~sknano.nanogen_gui._ng_model.NGModel`.

    """
    def __init__(self, args, model=None):
        app = QApplication(args)
        self.model = model
        self.view = NGView(self, self.model)
        self.model.notify_observers()
        self.view.show()
        app.exec_()

    def refresh_view(self):
        """Refresh `NGView`."""
        self.view.update_app_view()
