#! /usr/bin/env python
import sys
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import *

size = int(sys.argv[1])
N = size
mesh_size = size
source = (0, 0, 1, size)
sink = [(1, 0, size-1, 1), (size-1, 0, 1, size-1), (1, size-1, size-1, 1)]
conductance = 10**3  # this'll be used as resistance directly
hp = HeatProblem(N, source, sink, conductance, src_val=1.)

m = ROCModel(mesh_size)
m.load_problem(hp)
m.run_spice_solver()

print('Node Potentials')
print('---------------')
print_node_potentials(m)
print()
print('Node Currents')
print('-------------')
print_node_currents(m)
plot_heatmap(m, current_flow_plot='stream')
