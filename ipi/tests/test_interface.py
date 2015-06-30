"""Deals with testing the driver interface."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from ipi.interfaces.sockets import Client, Driver, InterfaceSocket


def test_client():
   """Client: startup without socket."""
   c = Client(_socket=False)


def test_driver():
   """Driver: startup without socket."""
   d = Driver(socket=None)


def test_interface():
   """InterfaceSocket: startup."""
   i = InterfaceSocket()
