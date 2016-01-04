"""Ensure that docs can always be built."""

# This file is part of i-PI.
# i-PI Copyright (C) 2015 i-PI developers
# See the "licenses" directory for full license information.

import os
import subprocess


def run_command(cmd):
    """Runs @cmd in doc directory."""
    cwd = os.getcwd()
    os.chdir(os.sep.join(__file__.split(os.sep)[:-1] + ["..", "doc"]))
    ret = subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    os.chdir(cwd)
    return ret


def test_make():
    """doc: run make"""
    ret = run_command("make")
    assert ret == 0


def test_make_aux():
    """doc: run make aux"""
    ret = run_command(["make", "aux"])
    assert ret == 0
