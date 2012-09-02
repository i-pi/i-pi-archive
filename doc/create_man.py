"""Help script which automatically generates the manual files.

Creates a latex file, corresponding to a section of the manual, for each of
the classes specified in help.py. It uses help.py to generate information 
about the tags for each class, and will include cross-references so that the
title of each tag corresponding to a different class will also be a hyperlink
in the manual to the section corresponding to that class.

Note that any new input class type must be added to the objects 
dictionary in help.py and the latex help file must be added to the end of 
the manual.lyx file for it to be included in the automatic help generation.

Also creates an xml file with the full list of all the tags.
"""

from help import help, objects
from help_list import help_list, list_objects

help(xml=True)
for opt in objects:
   help(latex=True, levels=1, option=opt, prefix=opt, ref=True)
for opt in list_objects:
   help_list(option=opt, prefix=opt, ref=True)
