#! /usr/bin/env python
import sys
import numpy as np
import itertools as it
from roc_model import ROCModel
import matplotlib.pyplot as plt
from heat_problem import HeatProblem


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

# save the input grid as well
# plt.imshow(grid, cmap='hot', interpolation='nearest')
# plt.savefig('tmp')


def print_current_table(m):
    frs = '{0:>9}    {1:>9}     {2}         {3:>1}'
    print(frs.format('Terminal1', 'Terminal2', 'Current', 'Direction'))
    checksum = 0.
    for mr in m.links:
        abs_current = abs(mr.ammeter.current)
        print(frs.format(str(mr.nodeblock1.coord),
                         str(mr.nodeblock2.coord),
                         '{:>10.6e}'.format(abs_current),
                         mr.cur_direction(mr.ammeter.current)))
        checksum += abs_current

    print('Checksum: {}'.format(checksum))

m.run_spice_solver(hp)

src_out = sum([n.sum_reduce_in_curs() for n in m.snk_nodeblocks()])
snk_in = sum([n.sum_reduce_out_curs() for n in m.src_nodeblocks()])

# for sanity checks
src_in = sum([n.sum_reduce_in_curs() for n in m.src_nodeblocks()])
snk_out = sum([n.sum_reduce_out_curs() for n in m.snk_nodeblocks()])
print('Source in : {}'.format(src_in))
print('Source out: {}'.format(src_out))
print('Sink in   : {}'.format(snk_in))
print('Sink out  : {}'.format(snk_out))
print_current_table(m)
fig, axes = plt.subplots(1,1)
potentials = np.array([[m.nodes[j][i].potential for i in
    range(mesh_size)] for j in range(mesh_size)])
axes.imshow(potentials, cmap='hot', interpolation='nearest')
U = np.array([[m.nodes[j][i].aggregate_current_vector()[0] for i in
        range(mesh_size)] for j in range(mesh_size)])
V = np.array([[m.nodes[j][i].aggregate_current_vector()[1] for i in
        range(mesh_size)] for j in range(mesh_size)])
# axes.streamplot(np.array([i for i in range(mesh_size)]),
               # np.array([i for i in range(mesh_size)]),
               # U, V, color='green', linewidth=1, density=0.5)

axes.quiver(np.array([i for i in range(mesh_size)]),
               np.array([i for i in range(mesh_size)]),
               U, V, color='green')

plt.show()
