#! /usr/bin/env python
import sys
import numpy as np
from roc_model import ROCModel
import matplotlib.pyplot as plt

# I am using python 3.6.1

# assume square grid
N = 100
num_iters = 1
pos = 50
mesh_size = int(sys.argv[1])
img_name = 'hm_{gr_sz}_{ms_sz}'

# initialize grid
grid = np.zeros((N, N))

# heat up an arbitrary point
grid[pos][pos] = 10.

cond_exp = -1
conductance = 10**cond_exp

sum_result = np.zeros((mesh_size,mesh_size))
# print(grid)
m = ROCModel(mesh_size)
m.load_problem(grid, conductance)

for i in range(num_iters):
    result = m.run_spice_solver()
    for s,r in zip(sum_result, result):
        s += r

debug = False
if debug:
    for i in range(N):
        for j in range(N):
            print('{0: >11.5f} '.format((result[i][j])*10.**15), end='')
        print()

# save the input grid as well
# plt.imshow(grid, cmap='hot', interpolation='nearest')
# plt.savefig('tmp')

def cross_plot(i, j):
    data_x = [sum_result[i][idx] for idx in range(mesh_size)]
    data_y = [sum_result[idx][j] for idx in range(mesh_size)]
    data_cross = [sum_result[idx][idx] for idx in range(mesh_size)]
    figs, axes = plt.subplots(1,3)
    axes[0].plot(range(mesh_size), data_x)
    axes[1].plot(range(mesh_size), data_y)
    axes[2].plot(range(mesh_size), data_cross)

def heat_plot():
    plt.imshow(sum_result, cmap='hot', interpolation='nearest')

tmp_x_loc = int(pos/int(N/mesh_size))
tmp_y_loc = int(pos/int(N/mesh_size))
cross_plot(tmp_x_loc,tmp_y_loc)
# plt.savefig(img_name.format(gr_sz=N, ms_sz=mesh_size))
plt.show()
