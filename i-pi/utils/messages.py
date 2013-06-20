"""Utility functions for outputting messages, diagnostics and errors'

Classes:
   Verbosity: Concise class to check the selected level of output

Functions:
   banner:    Prints the program welcome "screen"
   help:      Prints the input syntax help
   info:      Prints some information to standard output, depending on the level of verbosity
   warning:   Same as info, but with a "!W!" prefix and optionally printing a stack trace
"""

import traceback, sys

__all__ = ['Verbosity', 'verbosity',' help', 'banner', 'info', 'warning']


VERB_QUIET  = 0
VERB_LOW    = 1
VERB_MEDIUM = 2
VERB_HIGH   = 3
VERB_DEBUG  = 4

class Verbosity(object):

   level = "low"

   def __getattr__(self, name):
      if name is "quiet":
         return self.level >= VERB_QUIET
      elif name is "low":
         return self.level >= VERB_LOW
      elif name is "medium":
         return self.level >= VERB_MEDIUM
      elif name is "high":
         return self.level >= VERB_HIGH
      elif name is "debug":
         return self.level >= VERB_DEBUG

   def __setattr__(self, name, value):
      if name == "level":
         if value == "quiet":
            level = VERB_QUIET
         elif value == "low":
            level = VERB_LOW
         elif value == "medium":
            level = VERB_MEDIUM
         elif value == "high":
            level = VERB_HIGH
         elif value == "debug":
            level = VERB_DEBUG
         else: raise ValueError("Invalid verbosity level "+str(value)+" specified.")
         super(Verbosity,self).__setattr__("level", level)


verbosity = Verbosity()

def help():
   """Prints out a help string."""

   print """usage:  wrap-pi input """

def banner():
   """Prints out a banner."""

   print """
   ############################################################################
   #    _        _______  _____                                               #
   #   (_)      |_   __ \|_   _|    v. 0.9alpha                               #
   #   __  ______ | |__) | | |                                                #
   #  [  ||______||  ___/  | |      A python interface for (ab initio)        #
   #   | |       _| |_    _| |_     (path integral) molecular dynamics.       #
   #  [___]     |_____|  |_____|                                              #
   #                                                                          #
   ############################################################################
   """


def info(text="", show=True ):
   """Prints a warning message."""

   if not show:
      return
   print text

def warning(text="", show=True):
   """Prints a warning message."""
   if not show:
      return
   if verbosity.debug:
      traceback.print_stack(file=sys.stdout)
   print " !W! " + text
