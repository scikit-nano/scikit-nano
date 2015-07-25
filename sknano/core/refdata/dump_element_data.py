#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

from pkg_resources import resource_filename
import json
import os
import re
import sys

import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


_data_types = {'AtomicNumber': int, 'AtomicMass': float,
               'AtomicRadius': float, 'VanDerWaalsRadius': float,
               'SpaceGroupNumber': int}


def load_element_data():
    data = {}
    datafile = resource_filename('sknano', 'core/refdata/ElementData.dat')
    cre = re.compile(r'\d+\.\d*')
    with open(datafile, 'r') as f:
        headers = f.readline().strip().split('\t')
        for line in f:
            line = line.strip().split('\t')
            data[line[0]] = linedata = {}
            for k, v in zip(headers, line):
                if v == 'None':
                    linedata.update({k: None})
                    continue

                if k.startswith('Lattice'):
                    try:
                        v = [float(n) for n in cre.findall(v)]
                        linedata.update({k: v})
                    except ValueError as e:
                        print(e)
                        linedata.update({k: None})
                elif k in _data_types:
                    linedata.update({k: _data_types[k](v)})
                else:
                    linedata.update({k: v})
    return data


def dump_element_data(obj, fn, *args, **kwargs):
    with open(fn, 'wt') as fp:
        if fn.lower().endswith(("yaml", "yml")):
            if "Dumper" not in kwargs:
                kwargs["Dumper"] = Dumper
            yaml.dump(obj, fp, *args, **kwargs)
        else:
            json.dump(obj, fp, *args, **kwargs)


def main():
    element_data = load_element_data()
    dump_element_data(element_data,
                      os.path.join(os.pardir, 'element_data.yaml'))
    dump_element_data(element_data,
                      os.path.join(os.pardir, 'element_data.json'))


if __name__ == '__main__':
    sys.exit(main())
