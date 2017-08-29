"""Creates objects that deal with the alchemical exchanges."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
from ipi.utils.inputvalue import InputDictionary, InputAttribute, InputValue, InputArray, input_default


__all__ = ['InputAlchemy']


class InputAlchemy(InputDictionary):
    """Alchemy input class.

    Handles generating the appropriate alchemical exchange class from the xml input file.

    Fields:
        spicesA/B: isotopes for exchanges.
        timestep: An optional float giving the size of the timestep in atomic
            units. Defaults to 1.0.
    """
    attribs={"mode"  : (InputAttribute, {"dtype"   : str,
                                    "default" : 'dummy', 
                                    "help"    : " ",
                                    "options" : ['dummy']}) }

    fields = {
        "names"     : (InputArray, {"dtype"     : str,
                                    "default"   : input_default(factory=np.zeros, args=(0,), kwargs = {'dtype': np.dtype('|S6')}),
                                    "help"      : "The names of the atoms be to exchanged, in the format [name1, name2, ... ]." }),
        "nmc": (InputValue, {"dtype":     int,
                                  "default":   1,
                                  "help":      "The number of mc steps"})
             }

    dynamic = {  }

    default_help = "Holds all the information for doing Monte Carlo alchemical exchange moves. "
    default_label = "ALCHEMY"

    def store(self, alc):
        """Takes an alchemical exchange instance and stores a minimal representation of it.

        Args:
            alc: An alchemy object.
        """

        if alc == {}:
            return

        self.spicesA.store(alc.spicesA)
        self.spicesB.store(alc.spicesB)
        self.nmc.store(alc.nmc)

    def fetch(self):
        """Creates an ensemble object.

        Returns:
            An ensemble object of the appropriate mode and with the appropriate
            objects given the attributes of the InputEnsemble object.
        """

        rv = super(InputAlchemy, self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv
