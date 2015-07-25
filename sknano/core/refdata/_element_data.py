# -*- coding: utf-8 -*-
"""
================================================================
Reference data (:mod:`sknano.core.refdata._element_data`)
================================================================

.. currentmodule:: sknano.core.refdata._element_data

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

# from pkg_resources import resource_filename
import json
import os
# import re

import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


__all__ = ['element_data', '_load_element_data', '_dump_element_data']


def _load_element_data():
    data = None
    datafile = 'element_data.yaml'
    with open(os.path.join(os.path.dirname(__file__), datafile)) as fp:
        data = yaml.load(fp, Loader=Loader)
    return data


def _dump_element_data(obj, fn, *args, **kwargs):
    with open(fn, 'wt') as fp:
        if fn.lower().endswith(("yaml", "yml")):
            if "Dumper" not in kwargs:
                kwargs["Dumper"] = Dumper
            yaml.dump(obj, fp, *args, **kwargs)
        else:
            json.dump(obj, fp, *args, **kwargs)

element_data = _load_element_data()
