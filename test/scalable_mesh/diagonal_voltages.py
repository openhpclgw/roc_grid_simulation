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

exp_prob_size = 4
prob_size = 2**exp_prob_size
source = (0, 0, 1, 1)
sink = (prob_size-1, prob_size-1, 1, 1)
cond_exp = -3
conductance = 10**cond_exp
hp = HeatProblem(prob_size, source, sink, conductance, src_val=1.)

use_cached = True 

filename='tmp/diag_v_{}'

mesh = ROCModel(prob_size)

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

res_grid = get_results(mesh, cached=use_cached)

rect = 0.1,0.2,0.8,0.7
fig = plt.figure(figsize=(15,5))
ax = fig.add_axes(rect)

prob_range = range(prob_size)
data = [res_grid[i,i] for i in prob_range]

def custom_plot(ax, x, y):
    ax.plot(x, y,  marker='o',
                   markerfacecolor='none',
                   markeredgewidth=2)

    # put only 8 ticks on x axis
    ticks = [i for i in range(0, prob_size, int(prob_size/8))]

    ax.set_xlabel('Grid point')
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks)

    ax.set_ylabel('Potential')

    ax.set_ylim((0, max(data)*1.1))
    ax.set_xlim(0,prob_size-1)

    ax.grid(b=True, axis='x')
    ax.grid(b=True, axis='y', linestyle='dashed')

    ax.spines['top'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)

custom_plot(ax, prob_range, data)



plt.show()
