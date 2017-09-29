#! /usr/bin/env python
import sys
import numpy as np
from codegen import *
import matplotlib.pyplot as plt

# I am using python 3.6.1

# assume square grid
N = 100
mesh_size = int(sys.argv[1])
img_name = 'hm_{gr_sz}_{ms_sz}'

# initialize grid
grid = np.zeros((N,N))

# heat up an arbitrary point
grid[5][5] = 10.

cond_exp = -1
conductance = 10**cond_exp

# print(grid)
m = ROCModel(mesh_size)
m.load_problem(grid, conductance)
result = m.run_spice_solver()

debug = False
if debug:
    for i in range(N):
        for j in range(N):
            print('{0: >11.5f} '.format((result[i][j])*10.**15), end='')
        print()

# save the input grid as well
# plt.imshow(grid, cmap='hot', interpolation='nearest')
# plt.savefig('tmp')

plt.imshow(result, cmap='hot', interpolation='nearest')
plt.savefig(img_name.format(gr_sz=N, ms_sz=mesh_size))

