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

exp_largest_mesh = 4
prob_size = 2**exp_largest_mesh
source = (0, 0, 1, 1)
sink = (prob_size-1, prob_size-1, 1, 1)
cond_exp = -3
conductance = 10**cond_exp
hp = HeatProblem(prob_size, source, sink, conductance, src_val=10.)

use_cached = False

filename='tmp/bo_err_{}'

ground_truth_mesh = ROCModel(2**exp_largest_mesh)
exp_size_range = range(2,exp_largest_mesh)
meshes = {exp_size:ROCModel(2**exp_size) for exp_size in exp_size_range}


def blue_p(exp_mesh_size):
    val = 2**(exp_mesh_size-2)  # 1/4 point
    return (val,val)

def orange_p(exp_mesh_size):
    val = 3*2**(exp_mesh_size-2)  # 3/4 point
    return (val,val)

def get_results(mesh, cached=False):
    import os.path

    mesh.load_problem(hp)
    f = filename.format(mesh.w)
    exists = os.path.exists(f+'.out')
    if cached and exists:
        print('Using cache' + f)
        mesh.init_from_cache(f)
    else:
        mesh.run_spice_solver(f)
    return mesh.final_grid

base_grid = get_results(ground_truth_mesh, cached=use_cached)

blue_errs = []
orange_errs = []
for exp_mesh_size, mesh in meshes.items():
    print(exp_mesh_size)

    # grid is a numpy array
    tmp_grid = get_results(mesh, cached=use_cached)

    tmp_bp = blue_p(exp_mesh_size)
    tmp_op = orange_p(exp_mesh_size)
    blue_errs.append(abs(tmp_grid[tmp_bp]-base_grid[tmp_bp]))
    orange_errs.append(abs(tmp_grid[tmp_op]-base_grid[tmp_op]))

sizes = [2**i for i in exp_size_range]
print(blue_errs)
print(orange_errs)
fig, ax = plt.subplots(1,1)
ax.plot(sizes, blue_errs, sizes, orange_errs)
ax.set_xticks(sizes)
ax.set_xticklabels(sizes)
ax.set_ylim((0, max(max(blue_errs), max(orange_errs))*1.1))
ax.set_xlim(0,2**(exp_largest_mesh-1)*1.1)
plt.show()
