"""Ensure that docs can always be built."""

# This file is part of i-PI.
# i-PI Copyright (C) 2015 i-PI developers
# See the "licenses" directory for full license information.

import os
import subprocess


def test_make():
    """doc: run make"""
    cwd = os.getcwd()
    os.chdir(os.sep.join(__file__.split(os.sep)[:-1] + ["..", "doc"]))
    ret = subprocess.call("make")
    os.chdir(cwd)
    assert ret == 0
