#! /usr/bin/env python
import sys
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import *

size = int(sys.argv[1])
N = size
mesh_size = size
source = (0, 0, 1, size)
# sink = [(1, 0, size-1, 1), (size-1, 0, 1, size-1), (1, size-1, size-1, 1)]
sink = [(size-1, 0, 1, size-1)]
conductance = 10**-3  # this'll be used as resistance directly
hp = HeatProblem(N, source, sink, conductance, src_val=100.)

m = ROCModel(mesh_size)
m.load_problem(hp)
m.run_interconnect_solver(filename='interconnect_test')

# surface3d = True

# print('Node Potentials')
# print('---------------')
# print_node_potentials(m)
# print()
# print('Node Currents')
# print('-------------')
# print_node_currents(m)
# print()
# # print('Node Current Splits')
# # print('-------------')
# # print_current_splits(run_current_split_analysis(m))

# if not surface3d:
    # plot_heatmap(m, current_flow_plot=None)
# else:
    # plot_surface(m)
