"""Help script which automatically generates the manual files.

Also creates an xml file with the full list of all the tags.
"""

from help import help

help(xml=True)
help(latex=True, levels=1, option='simulation', prefix="simulation", ref=True)
help(latex=True, levels=1, option='atoms', prefix="atoms", ref=True)
help(latex=True, levels=1, option='barostats', prefix="barostats", ref=True)
help(latex=True, levels=1, option='beads', prefix="beads", ref=True)
help(latex=True, levels=1, option='cell', prefix="cell", ref=True)
help(latex=True, levels=1, option='forces', prefix="forces", ref=True)
help(latex=True, levels=1, option='ensembles', prefix="ensembles", ref=True)
help(latex=True, levels=1, option='interface', prefix="interface", ref=True)
help(latex=True, levels=1, option='thermostats', prefix="thermostats", ref=True)
help(latex=True, levels=1, option='prng', prefix="prng", ref=True)
