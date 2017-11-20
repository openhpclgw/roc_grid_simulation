#! /usr/bin/env python
import sys
import numpy as np
import itertools as it
import matplotlib as mpl
import matplotlib.pyplot as plt
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import (aggregate_current_vectors,
                            print_current_table,
                            energy_flow,
                            plot_heatmap)

exp_largest_mesh = 5
prob_size = 2**exp_largest_mesh
source = (0, 0, 1, 1)
sink = (prob_size-1, prob_size-1, 1, 1)
cond_exp = -3
conductance = 10**cond_exp
hp = HeatProblem(prob_size, source, sink, conductance, src_val=10.)

ground_truth_mesh = ROCModel(2**exp_largest_mesh)
exp_size_range = range(2,exp_largest_mesh)
meshes = {exp_size:ROCModel(2**exp_size) for exp_size in exp_size_range}

ground_truth_mesh.load_problem(hp)
ground_truth_mesh.run_spice_solver()
base_grid = ground_truth_mesh.final_grid

def blue_p(exp_mesh_size):
    val = 2**(exp_mesh_size-2)  # 1/4 point
    return (val,val)

def orange_p(exp_mesh_size):
    val = 3*2**(exp_mesh_size-2)  # 3/4 point
    return (val,val)


blue_errs = []
orange_errs = []
for exp_mesh_size, mesh in meshes.items():
    print(exp_mesh_size)
    mesh.load_problem(hp)
    mesh.run_spice_solver()

    # grid is a numpy array
    tmp_grid = mesh.final_grid

    tmp_bp = blue_p(exp_mesh_size)
    tmp_op = orange_p(exp_mesh_size)
    blue_errs.append(abs(tmp_grid[tmp_bp]-base_grid[tmp_bp]))
    orange_errs.append(abs(tmp_grid[tmp_op]-base_grid[tmp_op]))

sizes = [2**i for i in exp_size_range]
print(blue_errs)
print(orange_errs)
plt.plot(blue_errs, sizes, orange_errs, sizes)
plt.show()
