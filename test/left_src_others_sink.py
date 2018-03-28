#! /usr/bin/env python
import sys
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import *

def resistance(size, n1, n2):
    if n1[0]/size<0.5:
        return 0.002
    else:
        return 0.001

size = int(sys.argv[1])
N = size
mesh_size = size
source = (0, 0, 1, size)
# sink = [(2, 0, size-2, 1), (size-1, 0, 1, size-1), (2, size-1, size-2, 1)]
sink = [(size-1, 0, 1, size)]
hp = HeatProblem(N, source, sink, resistance)

norton = True

m = ROCModel(mesh_size, norton=norton)
m.load_problem(hp)
m.run_spice_solver()

surface3d = False

print('Node Potentials')
print('---------------')
print_node_potentials(m)
print()
print('Node Currents')
print('-------------')
print_node_currents(m)
print()
# print('Node Current Splits')
# print('-------------')
# splits = run_current_split_analysis(m)
# print_current_splits(splits)
# generate_sparams_from_splits(splits)

if norton:
    t = generate_loop_current_table(m)
    print(t)
    plot_heatmap_from_grid(t)
elif not surface3d:
    plot_heatmap(m, current_flow_plot=None)
else:
    plot_surface(m)
