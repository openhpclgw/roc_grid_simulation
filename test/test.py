#! /usr/bin/env python
import sys
import numpy as np
import itertools as it
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import (aggregate_current_vectors,
                            print_current_table,
                            energy_flow,
                            plot_surface,
                            plot_heatmap)


# I am using python 3.6.1

# assume square grid
# N = 100
# source = (10, 10, 5, 5)
# sink = (55, 10, 5, 50)
N = 20
source = (2, 2, 7, 10)
sink = [(11, 2, 1, 10), (5, 15, 2, 2)]
# source = (4, 4, 4, 12)
# sink = (12, 4, 4, 12)
cond_exp = -3
conductance = 10**cond_exp
mesh_size = int(sys.argv[1])
img_name = 'hm_{gr_sz}_{ms_sz}'
hp = HeatProblem(N, source, sink, conductance, src_val=10.)

surface3d = False  # if True 3d surface plot is generated

m = ROCModel(mesh_size)

m.load_problem(hp)
m.run_spice_solver()

eflow_data = energy_flow(m)
print('Source in : {}'.format(eflow_data['src_in']))
print('Source out: {}'.format(eflow_data['src_out']))
print('Sink in   : {}'.format(eflow_data['snk_in']))
print('Sink out  : {}'.format(eflow_data['snk_out']))

print_current_table(m)

# current flows currently doesn't work correctly with the new grid
# logic. It would probably require some post processing on the grid
# itself (which sucks)
if not surface3d:
    plot_heatmap(m)
else:
    plot_surface(m)
