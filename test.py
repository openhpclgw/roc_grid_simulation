#! /usr/bin/env python
import sys
import numpy as np
import itertools as it
from roc_model import ROCModel
import matplotlib.pyplot as plt

# I am using python 3.6.1

# assume square grid
N = 10
num_iters = 1
pos = 1
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
# m = ROCModel(mesh_size)
# m.load_problem(grid, conductance)

# for i in range(num_iters):
    # result = m.run_spice_solver()
    # for s,r in zip(sum_result, result):
        # s += r

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
    figs, axes = plt.subplots(3,3)
    axes[2][1].plot(range(mesh_size), data_x)
    axes[1][2].plot(range(mesh_size), data_y)
    axes[0][0].plot(range(mesh_size), data_cross)
    axes[1][1].imshow(sum_result, cmap='hot', interpolation='nearest')

def heat_plot(data):
    # zero out the source and the sink
    data[pos][pos] = 0.
    data[8][80] = 0.
    plt.imshow(data, cmap='hot', interpolation='nearest')

def numerical_solve(grid, sink, c, num_steps):
    grid2 = np.zeros((N,N))
    # for g,g2 in zip(grid,grid2):
    for (i,j) in it.product(range(N), range(N)):
        grid2[i][j] = grid[i][j]

    for step in range(num_steps):
        if step %2 == 0:
            for (i,j) in it.product(range(1, N-1), range(1, N-1)):
                if i==pos and j==pos:
                    continue
                if i==sink[0] and j==sink[1]:
                    continue
                grid2[i][j] = (grid[i+1][j] + grid[i-1][j] +
                                  grid[i][j+1] + grid[i][j-1] -
                                  0*grid[i][j])/c
            # print(grid2)
        else:
            for (i,j) in it.product(range(1, N-1), range(1, N-1)):
                if i==pos and j==pos:
                    continue
                if i==sink[0] and j==sink[1]:
                    continue
                grid[i][j] = (grid2[i+1][j] + grid2[i-1][j] +
                                  grid2[i][j+1] + grid2[i][j-1] -
                                  0*grid2[i][j])/c
            # print(grid)

        abs_delta = 0.
        for (i,j) in it.product(range(N), range(N)):
            abs_delta += abs(grid[i][j]-grid2[i][j])
        print(abs_delta)

    return grid

# tmp_x_loc = int(pos/int(N/mesh_size))
# tmp_y_loc = int(pos/int(N/mesh_size))
# cross_plot(tmp_x_loc,tmp_y_loc)
# plt.savefig(img_name.format(gr_sz=N, ms_sz=mesh_size))
heat_plot(numerical_solve(grid, (8,80), 4, 100))
plt.show()
