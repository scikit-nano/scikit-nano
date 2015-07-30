#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import os
import subprocess
import nose
from nose.tools import *


def test_gui_launch():
    if int(os.environ.get('SKNANO_SKIP_GUI_TESTS', 0)):
        print('Skipping nanogenui gui launch.')
    else:
        gui_launch_successful = False
        try:
            subprocess.call("nanogenui", timeout=1)
        except subprocess.TimeoutExpired:
            gui_launch_successful = True
        except OSError:
            pass
        assert_true(gui_launch_successful)


if __name__ == '__main__':
    nose.runmodule()
