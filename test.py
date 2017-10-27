#! /usr/bin/env python
import sys
import numpy as np
import itertools as it
from roc_model import ROCModel
import matplotlib.pyplot as plt
from heat_problem import HeatProblem
from analysis_utils import (aggregate_current_vectors,
                            print_current_table,
                            energy_flow,
                            node_potentials)


# I am using python 3.6.1

# assume square grid
# N = 100
# source = (10, 10, 5, 5)
# sink = (55, 10, 5, 50)
N = 20
source = (2, 2, 7, 10)
sink = (11, 2, 1, 10)
# source = (4, 4, 4, 12)
# sink = (12, 4, 4, 12)
cond_exp = -3
conductance = 10**cond_exp
num_iters = 1
pos = 1
mesh_size = int(sys.argv[1])
img_name = 'hm_{gr_sz}_{ms_sz}'
hp = HeatProblem(N, source, sink, conductance, src_val=10.)


m = ROCModel(mesh_size)

m.run_spice_solver(hp)

eflow_data = energy_flow(m)
print('Source in : {}'.format(eflow_data['src_in']))
print('Source out: {}'.format(eflow_data['src_out']))
print('Sink in   : {}'.format(eflow_data['snk_in']))
print('Sink out  : {}'.format(eflow_data['snk_out']))

print_current_table(m)

potentials = node_potentials(m)
fig, axes = plt.subplots(1,1)
axes.imshow(potentials, cmap='hot', interpolation='nearest')

U, V = aggregate_current_vectors(m)
# axes.streamplot(np.array([i for i in range(mesh_size)]),
               # np.array([i for i in range(mesh_size)]),
               # U, V, color='green', linewidth=1, density=0.5)

axes.quiver(np.array([i for i in range(mesh_size)]),
               np.array([i for i in range(mesh_size)]),
               U, V, color='green')

plt.show()
