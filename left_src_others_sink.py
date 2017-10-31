#! /usr/bin/env python
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import *

N = 5
mesh_size = 5
source = (0, 0, 1, 5)
sink = [(1, 0, 4, 1), (4, 0, 1, 4), (1, 4, 4, 1)]
cond_exp = -3
conductance = 10**cond_exp  # this'll be used as resistance directly
hp = HeatProblem(N, source, sink, conductance, src_val=10.)

m = ROCModel(mesh_size)
m.load_problem(hp)
m.run_spice_solver(hp)

print('Node Potentials')
print_node_potentials(m)
node_currents(m)
plot_heatmap(m, current_flow_plot='stream')
