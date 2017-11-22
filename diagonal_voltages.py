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

exp_prob_size = 5
prob_size = 2**exp_prob_size
source = (0, 0, 1, 1)
sink = (prob_size-1, prob_size-1, 1, 1)
cond_exp = -3
conductance = 10**cond_exp
hp = HeatProblem(prob_size, source, sink, conductance, src_val=10.)

use_cached = False

filename='tmp/diag_v_{}'

mesh = ROCModel(prob_size)

def get_results(mesh, cached=False):
    if cached:
        mesh.load_problem(hp)
        mesh.init_from_cache(filename.format(mesh.w))
    else:
        mesh.load_problem(hp)
        mesh.run_spice_solver(filename.format(mesh.w))
    return mesh.final_grid

res_grid = get_results(mesh, cached=use_cached)


fig, ax = plt.subplots(1,1)
prob_range = range(prob_size)
data = [res_grid[i,i] for i in prob_range]

ax.plot(prob_range, data)
ax.set_xticks(prob_range)
ax.set_xticklabels(prob_range)
ax.set_ylim((0, max(data)*1.1))
ax.set_xlim(0,prob_size-1)
plt.show()
