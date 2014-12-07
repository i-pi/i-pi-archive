"""Deals with testing the driver interface.
"""
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
