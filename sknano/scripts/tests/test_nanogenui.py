#! /usr/bin/env python

from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

import subprocess
import nose
from nose.tools import *
#from sknano.scripts.nanogenui import NanoGen


def test_gui_launch():
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
