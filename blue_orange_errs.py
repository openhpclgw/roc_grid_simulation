#! /usr/bin/env python
import sys
import numpy as np
import itertools as it
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import (aggregate_current_vectors,
                            print_current_table,
                            energy_flow,
                            node_potentials,
                            plot_heatmap)

exp_largest_mesh = 8
prob_size = 2**exp_largest_mesh
source = (0, 0, 1, 1)
sink = (prob_size, prob_size, 1, 1)
cond_exp = -3
conductance = 10**cond_exp
hp = HeatProblem(N, source, sink, conductance, src_val=10.)

ground_truth_mesh = ROCModel(2**exp_largest_mesh)
meshes = {exp_size:ROCModel(2**exp_size) for exp_size in
        range(2,exp_largest_mesh-1)}

ground_truth_mesh.load_problem(hp)
ground_truth_mesh.run_spice_solver()
base_grid = ground_truth_mesh.final_grid

def blue_p(exp_mesh_size):
    return 2**(exp_mesh_size-2)  # 1/4 point

def orange_p(exp_mesh_size):
    return 3*2**(exp_mesh_size-2)  # 3/4 point


blue_errs = []
orange_errs = []
for exp_mesh_size, mesh in meshes.items():
    mesh.load_problem(hp)
    mesh.run_spice_solver()
    tmp_grid = mesh.final_grid


eflow_data = energy_flow(m)
print('Source in : {}'.format(eflow_data['src_in']))
print('Source out: {}'.format(eflow_data['src_out']))
print('Sink in   : {}'.format(eflow_data['snk_in']))
print('Sink out  : {}'.format(eflow_data['snk_out']))

print_current_table(m)
plot_heatmap(m, current_flow_plot='stream')
