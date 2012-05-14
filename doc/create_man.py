import sys

src_dir = "../src/"

sys.path.append(src_dir)

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
