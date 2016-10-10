#!/usr/bin/env python2

import pytest
import mock
import tempfile

import ipi.engine.properties
import numpy as np



def create_a_fake_system_obj():

    natoms, q, atom_names, cell = 1, np.random.rand(3), ['Na'], np.random.rand(3,3)

    sys_attr = {
        'beads.natoms': natoms,
        'beads.q': q,
        'simul.step': 0,
        'beads.names': atom_names,
        'cell.h': cell
        }

    sys = mock.Mock(**sys_attr)
    # sys.beads.natoms = 'asdasd'

    stream = tempfile.NamedTemporaryFile(mode='w')

    return sys, stream

@mock.patch('ipi.engine.properties.io.print_file', autospec=True)
def test_Trajectories_print_traj(mock_io):

    system, stream = create_a_fake_system_obj()
    ipi.engine.properties.io.print_file = mock_io

    trj = ipi.engine.properties.Trajectories()
    trj.bind(system)


    trj.print_traj('positions', stream, b=0, format='xyz', cell_units='atomic_unit', flush=True )

    print mock_io.call_args_list
    #Just check that the arguments are the right ones!




test_Trajectories_print_traj()


# def test_create_a_fake_system_obj(monkeypatch):
#     natoms = 3
#     nbeads = 2

#     # per frames
#     q = ['random_number_'+str(x) for x in xrange(nbeads)]
#     qc = ['randoms_numbers']
#     f = ['random_number_'+str(x) for x in xrange(nbeads)]
#     temperature = 'random_number'

#     # need_to_sett = ['beads', 'cell']

#     # for _att in need_to_sett:
#     #     monkeypatch.setattr('ipi.engine.system', need_to_sett, raising=False)


#     set_values = {
#         'beads.natoms' : natoms,
#         'beads.nbeads' : nbeads,
#         'beads.q' : q,
#         'forces.f' : f,
#         'ensemble.temp': temperature,
#     }

#     for _key, _val in set_values.items():
#         attr_list = _key.split(('.'))
#         for _ii, _sk in enumerate(attr_list):
#             attr = '.'.join(attr_list[:_ii])

#             if not hasattr(ipi.engine.system, attr):
#                 if _ii == len(attr_list):
#                     monkeypatch.setattr('ipi.engine.system', attr, value=_val, raising=False)
#                 else:
#                     monkeypatch.setattr('ipi.engine.system', attr, raising=False)

#     print ipi.engine.system.beads.q
