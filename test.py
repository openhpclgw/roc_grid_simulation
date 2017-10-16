#! /usr/bin/env python
import sys
import numpy as np
from roc_model import ROCModel
import matplotlib.pyplot as plt

# I am using python 3.6.1

# assume square grid
N = 100
num_iters = 1
mesh_size = int(sys.argv[1])
img_name = 'hm_{gr_sz}_{ms_sz}'

# initialize grid
grid = np.zeros((N, N))

# heat up an arbitrary point
grid[5][5] = 10.

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

plt.imshow(sum_result, cmap='hot', interpolation='nearest')
# plt.savefig(img_name.format(gr_sz=N, ms_sz=mesh_size))
plt.show()
