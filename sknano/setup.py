#!/usr/bin/env python
from __future__ import division, print_function, absolute_import
from __future__ import unicode_literals


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('sknano', parent_package, top_path)
    config.add_subpackage('apps')
    config.add_subpackage('core')
    config.add_subpackage('generators')
    config.add_subpackage('io')
    config.add_subpackage('scripts')
    config.add_subpackage('structures')
    config.add_subpackage('testing')
    config.add_subpackage('utils')
    config.add_data_dir('data')
    #config.make_config_py()
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
